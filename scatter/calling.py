"""Library-aware, UMI-aware candidate calling for SCATTER.

The implementation deliberately keeps the evidence store in the requested
output directory.  It never creates an index beside, or otherwise writes to,
an input BAM.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import asdict
import csv
import math
from pathlib import Path
import sqlite3
from typing import Iterable

import pysam

from .annotation import load_genes as _load_genes, annotate_gene_overlap as _annotate_gene_overlap
from .clustering import SplitConfig, split_sites
from .io import read_barcodes
from .deeptssfilter import load_model, reverse_complement, score_sequences
from .protocol import (
    ASSAY_TYPES,
    LIBRARY_STRANDS,
    TSO_SOURCES,
    alignment_strand,
    biological_five_softclip,
    boolean_tag,
    clean_gene,
    evidence_five_position,
    evidence_five_softclip,
    has_five_softclip_alignment_anchor,
    load_gene_catalog,
    record_matches_mate,
    resolve_alignment_layout,
    resolve_evidence_mate,
    resolve_endpoint_rule,
    transcript_strand,
    transcript_five_position,
)


# This is the marker sequence used by CamoTSS and empirically observed in the
# K562 5' BAM.  It remains user-configurable because custom protocols can use
# another template-switch sequence.
DEFAULT_TSO_SEQUENCE = "TTTCTTATATGGG"
DEFAULT_MIN_GENE_READS = 50
DEFAULT_MIN_CLUSTER_READS = 50
DEFAULT_ADAPTIVE_GAP_MIN_BP = 15
DEFAULT_ADAPTIVE_GAP_MAX_BP = 100
DEFAULT_ADAPTIVE_GAP_PIVOT_READS = DEFAULT_MIN_CLUSTER_READS
DEFAULT_ADAPTIVE_GAP_DECAY = 0.30
DEFAULT_UTR5_FLANK_BP = 500
DEFAULT_OFF_KNOWN_TSS_BP = 100
DEFAULT_OFF_NEAR_TSS_BP = 500
DEFAULT_OFF_PROXIMAL_TSS_BP = 1000


def _load_barcode_whitelist(path: str | Path | None) -> set[str] | None:
    """Load the first non-empty TSV column, accepting a common header."""
    if path is None:
        return None
    return set(read_barcodes(path))


def _hamming_one_or_less(left: str, right: str) -> bool:
    return len(left) == len(right) and sum(a != b for a, b in zip(left, right)) <= 1


def _transcript_5_softclip(read: pysam.AlignedSegment, strand: str | None = None) -> str:
    """Compatibility wrapper for the selected read's physical 5′ soft clip."""
    return evidence_five_softclip(read)


def _scaled_mismatches(match_length: int, target_length: int, full_length_limit: int) -> int:
    """Scale a full-length mismatch allowance for a truncated observation."""
    return min(full_length_limit, math.floor(full_length_limit * match_length / target_length))


def _tso_match_length(
    read: pysam.AlignedSegment,
    tso_sequence: str,
    min_match: int,
    max_mismatches: int,
    *,
    technical_prefix_bp: int = 0,
    junction_slop_bp: int = 2,
    min_adjacent_match_bp: int = 6,
) -> int:
    """Return the supported TSO suffix length, or zero when no TSO is present.

    The known technical prefix is removed before matching. The remaining clip
    is normalized to sequencing orientation, and the TSO suffix must terminate
    at the alignment junction (allowing a small slop for untemplated bases).
    Full-length evidence is preferred; only genuinely shorter clips use the
    length-scaled partial-match rule.
    """
    if junction_slop_bp < 0:
        raise ValueError("junction_slop_bp must be non-negative")
    if not has_five_softclip_alignment_anchor(read, min_adjacent_match_bp):
        return 0
    clip = biological_five_softclip(read, technical_prefix_bp).upper()
    if read.is_reverse:
        clip = reverse_complement(clip)
    target = tso_sequence.upper()
    full_length_available = len(clip) >= len(target)
    for offset in range(junction_slop_bp + 1):
        available = len(clip) - offset
        # Do not rescue a failed full-length observation by silently dropping
        # junction bases and applying the more permissive truncated rule.
        if full_length_available and available < len(target):
            continue
        match_length = min(len(target), available)
        if match_length < min_match:
            continue
        observed = clip[available - match_length : available]
        expected = target[-match_length:]
        allowed = _scaled_mismatches(match_length, len(target), max_mismatches)
        if sum(left != right for left, right in zip(observed, expected)) <= allowed:
            return match_length
    return 0


def _has_tso(
    read: pysam.AlignedSegment,
    tso_sequence: str,
    min_match: int,
    max_mismatches: int,
    strand: str | None = None,
    *,
    technical_prefix_bp: int = 0,
    junction_slop_bp: int = 2,
    min_adjacent_match_bp: int = 6,
) -> bool:
    """Detect a junction-proximal TSO only in the biological 5′ soft clip."""
    return bool(_tso_match_length(
        read, tso_sequence, min_match, max_mismatches,
        technical_prefix_bp=technical_prefix_bp,
        junction_slop_bp=junction_slop_bp,
        min_adjacent_match_bp=min_adjacent_match_bp,
    ))


def _reference_matches_tso(
    fasta: pysam.FastaFile,
    chrom: str,
    read: pysam.AlignedSegment,
    tso_sequence: str,
    match_length: int,
    max_mismatches: int,
) -> bool:
    """Return whether adjacent genomic sequence can explain the TSO match.

    The queried side is determined by the selected record's physical 5′ end,
    not by the inferred transcript strand. This remains correct for reverse-
    stranded protocols and for explicit library-strand overrides.
    """
    if match_length < 1:
        return False
    try:
        if not read.is_reverse:
            end = int(read.reference_start)
            if end < match_length:
                return False
            sequence = fasta.fetch(chrom, end - match_length, end).upper()
        else:
            start = int(read.reference_end)
            sequence = reverse_complement(fasta.fetch(chrom, start, start + match_length).upper())
    except (KeyError, ValueError):
        return False
    if len(sequence) != match_length or "N" in sequence:
        return False
    expected = tso_sequence.upper()[-match_length:]
    allowed = _scaled_mismatches(match_length, len(tso_sequence), max_mismatches)
    return sum(left != right for left, right in zip(sequence, expected)) <= allowed


def _directional_umi_correction(umi_abundance: dict[str, int]) -> dict[str, str]:
    """Directional UMI correction (Hamming <=1) within one cell/cluster.

    A lower-abundance UMI is absorbed only by a strictly stronger neighbor,
    using the standard conservative ``low <= high / 2 - 1`` rule.  This avoids
    merging similarly abundant molecules merely because their UMIs are close.
    """
    ordered = sorted(umi_abundance, key=lambda umi: (-umi_abundance[umi], umi))
    parent = {umi: umi for umi in ordered}
    for high_index, high in enumerate(ordered):
        if parent[high] != high:
            continue
        high_count = umi_abundance[high]
        for low in ordered[high_index + 1 :]:
            if parent[low] != low:
                continue
            if umi_abundance[low] > (high_count / 2.0) - 1:
                continue
            if _hamming_one_or_less(high, low):
                parent[low] = high
    return parent


def _adaptive_gap(total_route_reads: int) -> int:
    """Return a read-depth-adaptive CTSS gap for one gene and route.

    Fifty route-qualified reads is the production cluster threshold and the
    pivot at which the full 100-bp gap is used.  More deeply observed genes
    receive progressively narrower gaps, bounded below by 15 bp.  Anchoring
    the curve at the production threshold makes the 100-bp upper bound
    effective for eligible low-depth genes rather than merely nominal.
    """
    relative_depth = max(1, total_route_reads) / DEFAULT_ADAPTIVE_GAP_PIVOT_READS
    gap = round(
        DEFAULT_ADAPTIVE_GAP_MAX_BP
        / (1 + DEFAULT_ADAPTIVE_GAP_DECAY * math.log10(relative_depth))
    )
    return max(DEFAULT_ADAPTIVE_GAP_MIN_BP, min(DEFAULT_ADAPTIVE_GAP_MAX_BP, gap))


def _weighted_median(counts: dict[int, int]) -> int:
    target = (sum(counts.values()) + 1) // 2
    running = 0
    for position in sorted(counts):
        running += counts[position]
        if running >= target:
            return position
    raise ValueError("Cannot calculate median of empty counts")


def _load_reference_tss(path: str | Path) -> dict[tuple[str, str, str], list[int]]:
    result: dict[tuple[str, str, str], list[int]] = defaultdict(list)
    with Path(path).open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result[(row["gene_id"].split(".")[0], row["chrom"], row["strand"])].append(int(row["tss_1"]))
    return {key: sorted(set(value)) for key, value in result.items()}


def _nearest(positions: list[int], value: int) -> int | None:
    if not positions:
        return None
    return min(abs(value - candidate) for candidate in positions)


def _write_tsv(path: Path, header: list[str], rows: Iterable[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


CLUSTER_COLUMNS = [
    "cluster_id", "gene_id", "chrom", "start_1", "end_1", "summit_1", "strand",
    "adaptive_gap_bp", "qualified_read_count", "raw_site_umi_count", "corrected_molecule_count", "ctss_count",
    "nearest_known_tss_distance_bp", "known_within_50bp", "deeptssfilter_score", "sequence_status", "model_sequence_length_bp",
    "tso_mode", "tso_matched_alignment_count", "tso_total_alignment_count", "tso_fraction_of_alignments",
    "candidate_source", "support_type", "tso_read_count", "off_read_count",
    "tso_cluster_count", "off_cluster_count", "source_cluster_count",
    "nearest_reference_tss_1", "signed_reference_tss_distance_bp",
    "off_tss_class", "off_tss_distance_pass", "off_tso_overlap",
    "off_route_filter_status", "off_route_filter_reason",
    "utr5_overlap", "utr5_flank_bp",
    "assay_type", "alignment_layout", "evidence_mate", "library_strand",
    "endpoint_rule", "tso_source", "tso_tag",
    "overlapping_gene_ids", "same_strand_overlapping_gene_ids", "opposite_strand_overlapping_gene_ids",
    "gene_assignment_status", "gene_assignment_evidence", "assigned_gene_span_overlap", "filter_status",
    "parent_cluster_id", "parent_start_1", "parent_end_1", "parent_qualified_read_count",
    "cluster_split_reason", "cluster_max_width_bp",
]
CTSS_COLUMNS = ["cluster_id", "gene_id", "chrom", "ctss_1", "strand", "corrected_molecule_count"]


def _load_utr5(path: str | Path, flank_bp: int) -> dict[tuple[str, str, str], list[tuple[int, int]]]:
    """Read, extend, and merge BED6 5′UTRs separately for each gene.

    Expanding before merging preserves the requested any-length overlap rule
    while turning large, redundant transcript-level annotations into a compact
    set of non-overlapping intervals for linear-time querying.
    """
    if flank_bp < 0:
        raise ValueError("utr5_flank_bp must be non-negative")
    result: dict[tuple[str, str, str], list[tuple[int, int]]] = defaultdict(list)
    with Path(path).open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6 or fields[5] not in {"+", "-"}:
                raise ValueError(f"5′UTR BED line {line_number} must have BED6 chrom/start/end/…/strand fields")
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"5′UTR BED line {line_number} has an invalid interval")
            gene = clean_gene(fields[3])
            # The v47 BED name is gene|transcript|5UTR. Accept a plain gene ID
            # as well, but never discard it: overlap must match the read's GX.
            if "|" in fields[3]:
                gene = clean_gene(fields[3].split("|", 1)[0])
            if gene is None:
                raise ValueError(f"5′UTR BED line {line_number} lacks a usable gene ID in column 4")
            result[(fields[0], fields[5], gene)].append((max(0, start - flank_bp), end + flank_bp))
    for key, intervals in result.items():
        intervals.sort()
        merged: list[tuple[int, int]] = []
        for start, end in intervals:
            if not merged or start > merged[-1][1]:
                merged.append((start, end))
            else:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        result[key] = merged
    if not result:
        raise ValueError("5′UTR BED contains no gene- and strand-resolved intervals")
    return result


def _mark_utr5_overlap(
    rows: list[dict[str, object]],
    utr5: dict[tuple[str, str, str], list[tuple[int, int]]],
    flank_bp: int,
) -> None:
    """Mark any-length overlap with the cluster's own GX gene 5′UTR."""
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["chrom"]), str(row["strand"]), str(row["gene_id"]))].append(row)
    for key, candidates in grouped.items():
        intervals = utr5.get(key, [])
        cursor = 0
        for row in sorted(candidates, key=lambda item: (int(item["start_1"]), int(item["end_1"]))):
            start_0, end_0 = int(row["start_1"]) - 1, int(row["end_1"])
            while cursor < len(intervals) and intervals[cursor][1] <= start_0:
                cursor += 1
            row["utr5_overlap"] = cursor < len(intervals) and intervals[cursor][0] < end_0
            row["utr5_flank_bp"] = flank_bp


def _annotate_reference_tss_support(
    rows: list[dict[str, object]],
    reference_tss: dict[tuple[str, str, str], list[int]],
    *,
    known_bp: int = DEFAULT_OFF_KNOWN_TSS_BP,
    near_bp: int = DEFAULT_OFF_NEAR_TSS_BP,
    proximal_bp: int = DEFAULT_OFF_PROXIMAL_TSS_BP,
) -> None:
    """Classify cluster summits by distance to the same gene's nearest TSS.

    Distance classes are inclusive at their upper boundaries.  The signed
    distance follows transcription direction: negative is upstream and
    positive is downstream on both genomic strands.  Retention itself uses
    the absolute distance, as requested for the three symmetric classes.
    """
    if not 0 <= known_bp <= near_bp <= proximal_bp:
        raise ValueError("off TSS distance thresholds must satisfy 0 <= known <= near <= proximal")
    for row in rows:
        gene, chrom, strand = str(row["gene_id"]), str(row["chrom"]), str(row["strand"])
        summit = int(row["summit_1"])
        positions = reference_tss.get((gene, chrom, strand), [])
        nearest = min(positions, key=lambda value: (abs(summit - value), value)) if positions else None
        if nearest is None:
            distance = signed = None
            category = "no_reference_tss"
            passed = False
        else:
            distance = abs(summit - nearest)
            signed = summit - nearest if strand == "+" else nearest - summit
            if distance <= known_bp:
                category = "known_tss"
                passed = True
            elif distance <= near_bp:
                category = "near_tss_shift"
                passed = True
            elif distance <= proximal_bp:
                category = "proximal_atss"
                passed = True
            else:
                category = "distal_filtered"
                passed = False
        row["nearest_reference_tss_1"] = nearest
        row["nearest_known_tss_distance_bp"] = distance
        row["signed_reference_tss_distance_bp"] = signed
        row["off_tss_class"] = category
        row["off_tss_distance_pass"] = passed


def _mark_off_tso_overlap(
    off_rows: list[dict[str, object]],
    tso_rows: list[dict[str, object]],
) -> None:
    """Mark inclusive interval overlap with a TSO cluster in the same locus."""
    tso_by_locus: dict[tuple[str, str, str], list[tuple[int, int]]] = defaultdict(list)
    for row in tso_rows:
        key = (str(row["chrom"]), str(row["gene_id"]), str(row["strand"]))
        tso_by_locus[key].append((int(row["start_1"]), int(row["end_1"])))
    for intervals in tso_by_locus.values():
        intervals.sort()

    grouped_off: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in off_rows:
        key = (str(row["chrom"]), str(row["gene_id"]), str(row["strand"]))
        grouped_off[key].append(row)
    for key, rows in grouped_off.items():
        intervals = tso_by_locus.get(key, [])
        cursor = 0
        for row in sorted(rows, key=lambda value: (int(value["start_1"]), int(value["end_1"]))):
            start, end = int(row["start_1"]), int(row["end_1"])
            while cursor < len(intervals) and intervals[cursor][1] < start:
                cursor += 1
            row["off_tso_overlap"] = bool(
                cursor < len(intervals) and intervals[cursor][0] <= end
            )


def _filter_off_candidates(
    off_rows: list[dict[str, object]],
    tso_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    """Retain reference-proximal off clusters or direct TSO overlaps."""
    _mark_off_tso_overlap(off_rows, tso_rows)
    retained: list[dict[str, object]] = []
    for row in off_rows:
        distance_pass = bool(row.get("off_tss_distance_pass"))
        tso_overlap = bool(row.get("off_tso_overlap"))
        if distance_pass and tso_overlap:
            row["off_route_filter_status"] = "retained"
            row["off_route_filter_reason"] = "reference_tss_window_and_tso_overlap"
        elif distance_pass:
            row["off_route_filter_status"] = "retained"
            row["off_route_filter_reason"] = str(row["off_tss_class"])
        elif tso_overlap:
            row["off_route_filter_status"] = "retained"
            row["off_route_filter_reason"] = "tso_overlap_rescue"
        else:
            row["off_route_filter_status"] = "filtered"
            row["off_route_filter_reason"] = str(row["off_tss_class"])
        if distance_pass or tso_overlap:
            retained.append(row)
    return retained


def _quantify_route_piece(conn, chrom, gene, strand, piece, gap, route):
    """Apply the unchanged per-cell UMI and CTSS rules to one interval."""
    start, end = piece[0][0], piece[-1][0]
    raw_count = sum(count for _, count in piece)
    rows = conn.execute(
        "SELECT cb,umi,pos,n_reads,tso_reads FROM evidence WHERE chrom=? AND gene=? AND strand=? AND pos BETWEEN ? AND ?",
        (chrom, gene, strand, start, end),
    ).fetchall()
    per_cell: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    per_position: dict[tuple[str, str], dict[int, int]] = defaultdict(lambda: defaultdict(int))
    total_alignments = tso_alignments = 0
    for cb, umi, pos, n_reads, tso_reads in rows:
        selected = tso_reads if route == "tso" else n_reads
        if not selected:
            continue
        per_cell[cb][umi] += selected
        per_position[(cb, umi)][pos] += selected
        total_alignments += n_reads
        tso_alignments += tso_reads
    corrected: Counter[int] = Counter()
    molecules = 0
    for cb, abundances in per_cell.items():
        parent = _directional_umi_correction(abundances)
        roots: dict[str, Counter[int]] = defaultdict(Counter)
        for umi, root in parent.items():
            roots[root].update(per_position[(cb, umi)])
        for counts in roots.values():
            corrected[min(counts, key=lambda pos: (-counts[pos], pos))] += 1
            molecules += 1
    if not corrected:
        return None
    return _route_piece_row(chrom, gene, strand, piece, gap, route, corrected, molecules,
                            total_alignments, tso_alignments)


def _route_piece_row(chrom, gene, strand, piece, gap, route, corrected, molecules,
                     total_alignments, tso_alignments):
    """One common schema for exact UMI counts and their proven-unique fast path."""
    start, end = piece[0][0], piece[-1][0]
    raw_count = sum(count for _, count in piece)
    selected_tso_reads = raw_count if route == "tso" else tso_alignments
    selected_off_reads = raw_count if route == "off" else 0
    return {
        "gene_id": gene, "chrom": chrom, "start_1": start, "end_1": end,
        "summit_1": _weighted_median(dict(corrected)), "strand": strand, "adaptive_gap_bp": gap,
        "qualified_read_count": raw_count, "raw_site_umi_count": raw_count,
        "corrected_molecule_count": molecules,
        "ctss_count": len(corrected), "tso_matched_alignment_count": tso_alignments,
        "tso_total_alignment_count": total_alignments,
        "tso_fraction_of_alignments": tso_alignments / total_alignments if total_alignments else 0.0,
        "candidate_source": route, "support_type": route,
        "tso_read_count": selected_tso_reads, "off_read_count": selected_off_reads,
        "tso_cluster_count": int(route == "tso"), "off_cluster_count": int(route == "off"),
        "source_cluster_count": 1,
        "utr5_overlap": None, "utr5_flank_bp": None, "_ctss": corrected,
    }


SPLIT_AUDIT_COLUMNS = [
    "parent_cluster_id", "route", "chrom", "gene_id", "strand", "start_1", "end_1",
    "parent_width_bp", "parent_reads", "retained_children", "discarded_children",
    "retained_reads", "discarded_reads", "split_reason", "max_child_width_bp",
]


def _subdivide_route_parent(conn, chrom, gene, strand, sites, gap, route,
                            parent_id, min_cluster_reads, config, split_audit):
    children = split_sites(sites, config, min_cluster_reads)
    parent_reads = sum(count for _, count in sites)
    retained = []
    discarded_reads = discarded_children = 0
    for piece, reason in children:
        reads = sum(count for _, count in piece)
        if reads < min_cluster_reads:
            discarded_reads += reads
            discarded_children += 1
            continue
        row = _quantify_route_piece(conn, chrom, gene, strand, piece, gap, route)
        if row is None:
            raise ValueError("A supported interval produced no UMI evidence")
        row.update(parent_cluster_id=parent_id, parent_start_1=sites[0][0],
                   parent_end_1=sites[-1][0], parent_qualified_read_count=parent_reads,
                   cluster_split_reason=reason, cluster_max_width_bp=config.max_width_bp)
        retained.append(row)
    kept_reads = sum(int(row["qualified_read_count"]) for row in retained)
    if kept_reads + discarded_reads != parent_reads:
        raise RuntimeError("Cluster subdivision failed read conservation")
    if split_audit is not None:
        split_audit.append(dict(
            parent_cluster_id=parent_id, route=route, chrom=chrom, gene_id=gene, strand=strand,
            start_1=sites[0][0], end_1=sites[-1][0], parent_width_bp=sites[-1][0]-sites[0][0]+1,
            parent_reads=parent_reads, retained_children=len(retained),
            discarded_children=discarded_children, retained_reads=kept_reads,
            discarded_reads=discarded_reads,
            split_reason="+".join(sorted({r for _, r in children})),
            max_child_width_bp=max((p[-1][0]-p[0][0]+1 for p, _ in children), default=0),
        ))
    return retained


def _route_candidates(
    conn: sqlite3.Connection,
    route: str,
    min_gene_reads: int = DEFAULT_MIN_GENE_READS,
    min_cluster_reads: int = DEFAULT_MIN_CLUSTER_READS,
    *, split_config: SplitConfig | None = None,
    split_audit: list[dict[str, object]] | None = None,
) -> list[dict[str, object]]:
    """Adaptive-gap components, valley/width subdivision, then original counting.

    Thresholds remain qualifying BAM-record counts. No reference or TSO
    rescue policy is changed by subdivision. Zero max width selects v0.4.4
    component boundaries for reproducibility.
    """
    if route not in {"off", "tso"}:
        raise ValueError(f"Unknown evidence route: {route}")
    if min_gene_reads < 1 or min_cluster_reads < 1:
        raise ValueError("Read thresholds must be positive integers")
    config = split_config or SplitConfig()
    expression = "tso_reads" if route == "tso" else "n_reads"
    candidates = []
    groups = conn.execute(
        f"SELECT chrom,gene,strand,SUM({expression}) FROM evidence "
        f"GROUP BY chrom,gene,strand HAVING SUM({expression}) >= ? ORDER BY chrom,gene,strand",
        (min_gene_reads,),
    ).fetchall()
    parent_number = 0
    for chrom, gene, strand, total_reads in groups:
        sites = conn.execute(
            f"SELECT pos,SUM({expression}) FROM evidence WHERE chrom=? AND gene=? AND strand=? "
            f"GROUP BY pos HAVING SUM({expression}) > 0 ORDER BY pos",
            (chrom, gene, strand),
        ).fetchall()
        gap = _adaptive_gap(int(total_reads))
        pieces = []
        for site in sites:
            if not pieces or site[0] - pieces[-1][-1][0] > gap:
                pieces.append([site])
            else:
                pieces[-1].append(site)
        for sites in pieces:
            if sum(count for _, count in sites) < min_cluster_reads:
                continue
            parent_number += 1
            candidates.extend(_subdivide_route_parent(
                conn, chrom, gene, strand, sites, gap, route,
                f"{route.upper()}_PARENT_{parent_number:08d}", min_cluster_reads, config, split_audit,
            ))
    return candidates

def _merge_mix_candidates(
    tso: list[dict[str, object]],
    off: list[dict[str, object]],
) -> list[dict[str, object]]:
    """Union overlapping TSO and retained all-read candidates.

    The all-read (``off``) route already contains its overlapping TSO reads, so
    its representative is retained whenever a component contains one. This
    avoids double-counting evidence. Components without an all-read candidate
    retain their strongest TSO candidate. Source annotations distinguish
    ``tso``, ``off``, and ``mix_tso+off`` support.
    """
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in tso + off:
        grouped[(str(row["chrom"]), str(row["gene_id"]), str(row["strand"]))].append(row)
    result: list[dict[str, object]] = []
    for rows in grouped.values():
        pieces: list[list[dict[str, object]]] = []
        for row in sorted(rows, key=lambda item: (int(item["start_1"]), int(item["end_1"]))):
            if not pieces or int(row["start_1"]) > max(int(item["end_1"]) for item in pieces[-1]):
                pieces.append([row])
            else:
                pieces[-1].append(row)
        for piece in pieces:
            tso_members = [row for row in piece if row["candidate_source"] == "tso"]
            off_members = [row for row in piece if row["candidate_source"] == "off"]
            chosen = max(off_members or piece, key=lambda row: int(row["corrected_molecule_count"]))
            merged = dict(chosen)
            source_set = {str(row["candidate_source"]) for row in piece}
            if source_set == {"tso", "off"}:
                merged["candidate_source"] = "mix_tso+off"
                merged["support_type"] = "both"
            else:
                merged["candidate_source"] = next(iter(source_set))
                merged["support_type"] = merged["candidate_source"]
            merged["tso_cluster_count"] = len(tso_members)
            merged["off_cluster_count"] = len(off_members)
            merged["source_cluster_count"] = len(piece)
            result.append(merged)
    return result


def _add_protocol_metadata(
    rows: Iterable[dict[str, object]],
    *,
    assay_type: str,
    alignment_layout: str,
    evidence_mate: str,
    library_strand: str,
    tso_source: str,
    tso_tag: str,
    endpoint_rule: str = "physical-5prime",
) -> None:
    position_label = ("selected_alignment_transcript_5prime" if endpoint_rule == "transcript-5prime"
                      else "selected_read_physical_5prime")
    for row in rows:
        row["assay_type"] = assay_type
        row["alignment_layout"] = alignment_layout
        row["evidence_mate"] = evidence_mate
        row["library_strand"] = library_strand
        row["endpoint_rule"] = position_label
        row["tso_source"] = tso_source
        row["tso_tag"] = tso_tag if tso_source == "tag" else ""


def _assign_route_cluster_ids(
    rows: Iterable[dict[str, object]], prefix: str
) -> None:
    """Assign deterministic IDs before route rows are copied or exported."""
    for number, row in enumerate(rows, start=1):
        row["cluster_id"] = f"{prefix}_{number:08d}"


def _output_paths(root: Path, tso_mode: str) -> dict[str, Path]:
    """Return explicit, route-aware paths for one SCATTER call."""
    clusters_dir = root / "clusters"
    paths = {
        "final_raw_clusters": clusters_dir / f"{tso_mode}_route_final_raw_clusters.tsv",
    }
    if tso_mode in {"tso", "mix"}:
        paths["tso_route_raw_clusters"] = clusters_dir / "tso_route_raw_clusters.tsv"
    if tso_mode in {"off", "mix"}:
        paths.update({
            "off_route_before_reference_tss_filter_raw_clusters": (
                clusters_dir / "off_route_before_reference_tss_filter_raw_clusters.tsv"
            ),
            "off_route_after_reference_tss_filter_raw_clusters": (
                clusters_dir / "off_route_after_reference_tss_filter_raw_clusters.tsv"
            ),
        })
    if tso_mode == "mix":
        paths["mix_route_off_after_tso_rescue_raw_clusters"] = (
            clusters_dir / "mix_route_off_after_tso_rescue_raw_clusters.tsv"
        )
    return paths


def call_clusters(
    *, bam_path: str,
    bam_index: str,
    fasta_path: str,
    reference_tss_path: str,
    model_path: str,
    outdir: str,
    library_type: str | None = None,
    alignment_layout: str | None = None,
    assay_type: str = "5prime",
    evidence_mate: str | None = "auto",
    library_strand: str = "gene",
    endpoint_rule: str = "auto",
    tso_source: str = "evidence",
    tso_tag: str = "ZT",
    cb_tag: str = "CB",
    umi_tag: str = "UB",
    barcode_whitelist: str | None = None,
    min_mapq: int = 30,
    min_gene_reads: int = DEFAULT_MIN_GENE_READS,
    min_cluster_reads: int = DEFAULT_MIN_CLUSTER_READS,
    threshold: float = 0.9,
    tso_mode: str = "off",
    tso_sequence: str = DEFAULT_TSO_SEQUENCE,
    tso_min_match: int = 8,
    tso_max_mismatches: int = 2,
    evidence_technical_prefix_bp: int = 0,
    tso_junction_slop_bp: int = 2,
    tso_min_adjacent_match_bp: int = 6,
    tso_reject_genomic_match: bool = True,
    utr5_bed: str | None = None,
    utr5_flank_bp: int = DEFAULT_UTR5_FLANK_BP,
    off_known_tss_bp: int = DEFAULT_OFF_KNOWN_TSS_BP,
    off_near_tss_bp: int = DEFAULT_OFF_NEAR_TSS_BP,
    off_proximal_tss_bp: int = DEFAULT_OFF_PROXIMAL_TSS_BP,
    contigs: set[str] | None = None,
    max_records: int | None = None,
    raw_only: bool = False,
    split_config: SplitConfig | None = None,
    device: str = "auto",
    batch_size: int = 256,
    threads: int = 4,
    sample_id: str | None = None,
) -> dict[str, object]:
    """Call raw clusters and DeepTSSFilter-passing clusters.

    Layout, evidence mate, transcript strand, and TSO source are independent.
    Paired input can select Read 1 or Read 2. ``library_strand=forward`` keeps
    alignment strand, ``reverse`` flips it, and ``gene`` uses the unique GX
    gene's annotated strand. Single-read 5-prime input defaults to the
    alignment boundary facing transcript 5-prime: left on RNA '+' and right
    on RNA '-'. Other layouts retain the physical read 5-prime boundary.
    Explicit endpoint_rule='physical-5prime' reproduces the earlier single-R2
    convention. TSO soft clips remain tied to the physical sequenced end.
    A transferred boolean BAM tag can supply TSO status
    when the TSO-bearing mate is absent from the coordinate evidence record.

    ``off`` clusters all commonly qualified reads without performing TSO
    validation, then retains off-only summits no farther than 1 kb from the
    same gene's nearest reference TSS. ``tso`` clusters only TSO-positive
    reads without a reference-distance restriction. ``mix`` also rescues
    any off cluster overlapping a TSO cluster before merging route candidates.
    """
    split_config = split_config or SplitConfig()
    if not 0 <= threshold <= 1:
        raise ValueError('Model threshold must be between 0 and 1')
    if tso_mode == "on":  # Backward-compatible spelling for the TSO-only route.
        tso_mode = "tso"
    resolved_layout = resolve_alignment_layout(library_type, alignment_layout)
    resolved_mate = resolve_evidence_mate(resolved_layout, evidence_mate)
    resolved_endpoint_rule = resolve_endpoint_rule(endpoint_rule, assay_type, resolved_layout)
    if assay_type not in ASSAY_TYPES:
        raise ValueError(f"Unknown assay type: {assay_type}")
    if library_strand not in LIBRARY_STRANDS:
        raise ValueError(f"Unknown library strand: {library_strand}")
    if tso_source not in TSO_SOURCES:
        raise ValueError(f"Unknown TSO source: {tso_source}")
    if tso_mode == "mixed":
        tso_mode = "mix"  # Backward compatibility for direct Python callers.
    if tso_mode not in {"off", "tso", "mix"}:
        raise ValueError(f"Unknown TSO mode: {tso_mode}")
    if tso_mode in {"tso", "mix"} and tso_source == "none":
        raise ValueError("TSO mode 'tso' or 'mix' requires --tso-source evidence or tag")
    effective_tso_source = tso_source if tso_mode in {"tso", "mix"} else "none"
    if tso_source == "tag" and len(tso_tag) != 2:
        raise ValueError("BAM TSO tag must be exactly two characters")
    if utr5_flank_bp < 0:
        raise ValueError("utr5_flank_bp must be non-negative")
    if not 0 <= off_known_tss_bp <= off_near_tss_bp <= off_proximal_tss_bp:
        raise ValueError("off TSS distance thresholds must satisfy 0 <= known <= near <= proximal")
    if min_gene_reads < 1 or min_cluster_reads < 1:
        raise ValueError("Read thresholds must be positive integers")
    tso_sequence = tso_sequence.upper()
    if not tso_sequence or set(tso_sequence).difference("ACGT"):
        raise ValueError("TSO sequence must contain only A/C/G/T")
    if not 1 <= tso_min_match <= len(tso_sequence):
        raise ValueError("TSO minimum match must be between 1 and the TSO sequence length")
    if tso_max_mismatches < 0 or tso_max_mismatches >= tso_min_match:
        raise ValueError("TSO maximum mismatches must be non-negative and below the minimum match length")
    if evidence_technical_prefix_bp < 0:
        raise ValueError("Evidence technical prefix length must be non-negative")
    if tso_junction_slop_bp < 0:
        raise ValueError("TSO junction slop must be non-negative")
    if tso_min_adjacent_match_bp < 1:
        raise ValueError("TSO adjacent alignment anchor must be positive")
    gene_catalog = load_gene_catalog(reference_tss_path)
    allowed_barcodes = _load_barcode_whitelist(barcode_whitelist)
    root = Path(outdir).expanduser().resolve()
    clusters_dir = root / "clusters"
    ctss_dir = root / "ctss"
    work_dir = root / "work"
    output_paths = _output_paths(root, tso_mode)
    output_directories = [clusters_dir, work_dir]
    if not raw_only:
        output_directories.append(ctss_dir)
    for directory in output_directories:
        directory.mkdir(parents=True, exist_ok=True)
    db_path = work_dir / "ctss_evidence.sqlite"
    if db_path.exists():
        raise FileExistsError(f"Refusing to overwrite an existing evidence database: {db_path}")

    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute(
        "CREATE TABLE evidence (chrom TEXT, gene TEXT, strand TEXT, pos INTEGER, cb TEXT, umi TEXT, "
        "n_reads INTEGER NOT NULL, tso_reads INTEGER NOT NULL, "
        "PRIMARY KEY(chrom, gene, strand, pos, cb, umi)) WITHOUT ROWID"
    )
    sql = (
        "INSERT INTO evidence VALUES (?, ?, ?, ?, ?, ?, 1, ?) "
        "ON CONFLICT(chrom,gene,strand,pos,cb,umi) DO UPDATE SET "
        "n_reads=n_reads+1, tso_reads=tso_reads+excluded.tso_reads"
    )
    counters = Counter()
    pending: list[tuple[object, ...]] = []
    genomic_tso_cache: dict[tuple[str, bool, int, int, int], bool] = {}
    genomic_tso_cache_limit = 1_000_000
    tso_fasta = (
        pysam.FastaFile(fasta_path)
        if tso_mode in {"tso", "mix"} and tso_source == "evidence" and tso_reject_genomic_match
        else None
    )
    with pysam.AlignmentFile(bam_path, "rb", index_filename=bam_index) as bam:
        if contigs is None:
            reads = bam.fetch(until_eof=True)
        else:
            unavailable = sorted(contigs.difference(bam.references))
            if unavailable:
                raise ValueError(f"Requested contigs are absent from BAM: {', '.join(unavailable)}")
            reads = (read for chrom in sorted(contigs) for read in bam.fetch(chrom))
        for read in reads:
            if max_records is not None and counters["records_examined"] >= max_records:
                break
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            counters["records_examined"] += 1
            if not record_matches_mate(read, resolved_mate):
                continue
            counters["selected_mate_primary_alignments"] += 1
            if read.mapping_quality < min_mapq:
                counters["low_mapq"] += 1
                continue
            if not (read.has_tag(cb_tag) and read.has_tag(umi_tag) and read.has_tag("GX")):
                counters["missing_required_tag"] += 1
                continue
            gene = clean_gene(read.get_tag("GX"))
            if gene is None:
                counters["ambiguous_or_unassigned_gene"] += 1
                continue
            gene_record = gene_catalog.get(gene)
            if gene_record is None:
                counters["gene_absent_from_reference"] += 1
                continue
            chrom = bam.get_reference_name(read.reference_id)
            if chrom != gene_record.chrom:
                counters["gene_chromosome_mismatch"] += 1
                continue
            cb, umi = read.get_tag(cb_tag), read.get_tag(umi_tag)
            if not cb or not umi or cb == "-" or umi == "-":
                counters["invalid_cb_or_umi"] += 1
                continue
            if allowed_barcodes is not None and cb not in allowed_barcodes:
                counters["barcode_not_whitelisted"] += 1
                continue
            counters["candidate_alignments_before_tso_filter"] += 1
            observed_strand = alignment_strand(read)
            strand = transcript_strand(read, gene_record, library_strand)
            counters["alignment_same_as_gene_strand"] += int(observed_strand == gene_record.strand)
            counters["alignment_opposite_gene_strand"] += int(observed_strand != gene_record.strand)
            counters["resolved_strand_matches_gene"] += int(strand == gene_record.strand)
            # OFF is deliberately TSO-agnostic. TSO status is calculated only
            # for the modes that split or validate evidence by TSO status.
            if tso_mode not in {"tso", "mix"}:
                tso_match = 0
            elif tso_source == "evidence":
                match_length = _tso_match_length(
                    read, tso_sequence, tso_min_match, tso_max_mismatches,
                    technical_prefix_bp=evidence_technical_prefix_bp,
                    junction_slop_bp=tso_junction_slop_bp,
                    min_adjacent_match_bp=tso_min_adjacent_match_bp,
                )
                counters["tso_sequence_matched_alignments"] += int(bool(match_length))
                genomic_match = False
                if match_length and tso_fasta is not None:
                    cache_key = (
                        chrom, bool(read.is_reverse), int(read.reference_start),
                        int(read.reference_end), match_length,
                    )
                    cached_match = genomic_tso_cache.get(cache_key)
                    if cached_match is None:
                        if len(genomic_tso_cache) >= genomic_tso_cache_limit:
                            genomic_tso_cache.clear()
                            counters["tso_genomic_cache_resets"] += 1
                        cached_match = _reference_matches_tso(
                            tso_fasta, chrom, read, tso_sequence,
                            match_length, tso_max_mismatches,
                        )
                        genomic_tso_cache[cache_key] = cached_match
                        counters["tso_genomic_cache_misses"] += 1
                    else:
                        counters["tso_genomic_cache_hits"] += 1
                    genomic_match = cached_match
                counters["tso_genomic_match_rejected_alignments"] += int(genomic_match)
                tso_match = int(bool(match_length) and not genomic_match)
            else:
                try:
                    tso_match = int(boolean_tag(read, tso_tag))
                except KeyError as exc:
                    raise ValueError(
                        f"Selected evidence record lacks required transferred TSO tag {tso_tag}: {read.query_name}"
                    ) from exc
            position = (transcript_five_position(read, strand) if resolved_endpoint_rule == "transcript-5prime"
                        else evidence_five_position(read))
            pending.append((chrom, gene, strand, position, cb, umi, tso_match))
            counters["eligible_alignments"] += 1
            counters["tso_matched_alignments"] += tso_match
            if len(pending) >= 50_000:
                conn.executemany(sql, pending)
                conn.commit()
                pending.clear()
    if tso_fasta is not None:
        tso_fasta.close()
    if pending:
        conn.executemany(sql, pending)
        conn.commit()
    if not counters["selected_mate_primary_alignments"]:
        conn.close()
        raise ValueError(
            f"Selected evidence mate {resolved_mate} has no primary mapped alignments; "
            "choose the BAM mate carrying the cDNA alignment"
        )
    if not counters["eligible_alignments"]:
        conn.close()
        raise ValueError(
            f"Selected evidence mate {resolved_mate} has no alignments passing MAPQ/GX/CB/UB/reference checks"
        )
    if library_strand != "gene" and counters["candidate_alignments_before_tso_filter"] >= 1000:
        strand_agreement = counters["resolved_strand_matches_gene"] / counters["candidate_alignments_before_tso_filter"]
        if strand_agreement < 0.95:
            conn.close()
            raise ValueError(
                f"Configured --library-strand {library_strand} agrees with GX gene strand for only "
                f"{strand_agreement:.1%} of qualified alignments; expected at least 95%"
            )
    counters["unique_site_cell_umi_observations"] = conn.execute("SELECT COUNT(*) FROM evidence").fetchone()[0]

    ref_tss = _load_reference_tss(reference_tss_path)
    split_audit: list[dict[str, object]] = []
    tso_candidates = _route_candidates(
        conn, "tso", min_gene_reads, min_cluster_reads,
        split_config=split_config, split_audit=split_audit,
    ) if tso_mode in {"tso", "mix"} else []
    off_before_tss_filter = _route_candidates(
        conn, "off", min_gene_reads, min_cluster_reads,
        split_config=split_config, split_audit=split_audit,
    ) if tso_mode in {"off", "mix"} else []
    _annotate_reference_tss_support(
        tso_candidates + off_before_tss_filter,
        ref_tss,
        known_bp=off_known_tss_bp,
        near_bp=off_near_tss_bp,
        proximal_bp=off_proximal_tss_bp,
    )
    # Assign stable route-level IDs before making the standalone-off copies.
    # These IDs make every exported raw-route table independently auditable;
    # the final mix table receives its own SCATTER_* IDs below.
    _assign_route_cluster_ids(tso_candidates, "TSO_RAW")
    _assign_route_cluster_ids(off_before_tss_filter, "OFF_PRETSS")
    # Keep the standalone off/notTSO route independent of TSO evidence.  A
    # separate copy is allowed to use TSO-overlap rescue only inside mix
    # mode; otherwise exporting notTSO from a mix scan would silently make
    # that route TSO-dependent.
    off_reference_candidates = _filter_off_candidates(
        [dict(row) for row in off_before_tss_filter], []
    )
    off_mix_candidates = (
        _filter_off_candidates(off_before_tss_filter, tso_candidates)
        if tso_mode == "mix"
        else off_reference_candidates
    )
    if tso_mode == "tso":
        raw_clusters = tso_candidates
    elif tso_mode == "off":
        raw_clusters = off_reference_candidates
    else:
        raw_clusters = _merge_mix_candidates(tso_candidates, off_mix_candidates)
    _add_protocol_metadata(
        list({
            id(row): row for row in (
                tso_candidates + off_before_tss_filter + off_reference_candidates
                + off_mix_candidates + raw_clusters
            )
        }.values()),
        assay_type=assay_type,
        alignment_layout=resolved_layout,
        evidence_mate=resolved_mate,
        library_strand=library_strand,
        tso_source=effective_tso_source,
        tso_tag=tso_tag,
        endpoint_rule=resolved_endpoint_rule,
    )
    if tso_mode in {"tso", "mix"}:
        counters["tso_route_raw_clusters"] = len(tso_candidates)
    if tso_mode in {"off", "mix"}:
        counters["off_route_raw_clusters_before_tss_filter"] = len(off_before_tss_filter)
        counters["off_route_clusters_retained_by_reference_tss"] = len(off_reference_candidates)
        counters["off_route_clusters_filtered_by_tss_distance"] = (
            len(off_before_tss_filter) - len(off_reference_candidates)
        )
    if tso_mode == "mix":
        counters["mix_off_clusters_after_tso_rescue"] = len(off_mix_candidates)
        counters["mix_off_clusters_rescued_by_tso_overlap"] = sum(
            row.get("off_route_filter_reason") == "tso_overlap_rescue"
            for row in off_mix_candidates
        )

    gene_intervals = _load_genes(reference_tss_path)
    _annotate_gene_overlap(
        list({id(row): row for row in (
            tso_candidates + off_before_tss_filter + off_reference_candidates
            + off_mix_candidates + raw_clusters
        )}.values()), gene_intervals,
    )
    audit_path = root / "metadata" / "cluster_splitting_audit.tsv"
    audit_path.parent.mkdir(exist_ok=True)
    _write_tsv(audit_path, SPLIT_AUDIT_COLUMNS, split_audit)
    output_paths["cluster_splitting_audit"] = audit_path

    # Keep the all-read candidates before and after the reference-TSS selection
    # auditable, without passing their pre-selection rows to DeepTSSFilter.
    for row in off_before_tss_filter:
        row["tso_mode"] = tso_mode
        row.setdefault("filter_status", "pending_reference_tss_filter")
    if tso_mode in {"off", "mix"}:
        _write_tsv(
            output_paths["off_route_before_reference_tss_filter_raw_clusters"],
            CLUSTER_COLUMNS,
            off_before_tss_filter,
        )
        _write_tsv(
            output_paths["off_route_after_reference_tss_filter_raw_clusters"],
            CLUSTER_COLUMNS,
            off_reference_candidates,
        )
    if tso_mode == "mix":
        _write_tsv(
            output_paths["mix_route_off_after_tso_rescue_raw_clusters"],
            CLUSTER_COLUMNS,
            off_mix_candidates,
        )
    if tso_mode in {"tso", "mix"}:
        _write_tsv(
            output_paths["tso_route_raw_clusters"],
            CLUSTER_COLUMNS,
            tso_candidates,
        )

    if raw_only:
        for cluster_number, row in enumerate(raw_clusters, start=1):
            row["cluster_id"] = f"SCATTER_{cluster_number:08d}"
            row["tso_mode"] = tso_mode
            row["filter_status"] = "raw_only"
        conn.close()
        _write_tsv(output_paths["final_raw_clusters"], CLUSTER_COLUMNS, raw_clusters)
        return {
            "cluster_splitting": asdict(split_config),
            "evidence_database": str(db_path), "raw_clusters": len(raw_clusters),
            "raw_only": True, "tso_mode": tso_mode,
            "endpoint_rule": resolved_endpoint_rule, "requested_endpoint_rule": endpoint_rule,
            "output_files": {key: str(path) for key, path in output_paths.items()},
            "tso_sequence": tso_sequence if tso_mode in {"tso", "mix"} and tso_source == "evidence" else None,
            "assay_type": assay_type, "alignment_layout": resolved_layout,
            "evidence_mate": resolved_mate, "library_strand": library_strand,
            "tso_source": effective_tso_source,
            "configured_tso_source": tso_source,
            "tso_tag": tso_tag if effective_tso_source == "tag" else None,
            "evidence_technical_prefix_bp": evidence_technical_prefix_bp,
            "tso_junction_slop_bp": tso_junction_slop_bp,
            "tso_min_adjacent_match_bp": tso_min_adjacent_match_bp,
            "tso_reject_genomic_match": tso_reject_genomic_match,
            "min_gene_reads": min_gene_reads, "min_cluster_reads": min_cluster_reads,
            "off_known_tss_bp": off_known_tss_bp if tso_mode in {"off", "mix"} else None,
            "off_near_tss_bp": off_near_tss_bp if tso_mode in {"off", "mix"} else None,
            "off_proximal_tss_bp": off_proximal_tss_bp if tso_mode in {"off", "mix"} else None,
            **dict(counters),
        }

    from .model import score_loci
    scores, sequence_statuses, inference_metadata = score_loci(
        raw_clusters, fasta_path, model_path, device=device, batch_size=batch_size, threads=threads,
    )
    ctss_rows: list[dict[str, object]] = []
    for cluster_number, row in enumerate(raw_clusters, start=1):
        row["cluster_id"] = f"SCATTER_{cluster_number:08d}"
        row["tso_mode"] = tso_mode
        summit, chrom, strand = int(row["summit_1"]), str(row["chrom"]), str(row["strand"])
        locus = (chrom, summit, strand)
        score = scores[locus]
        row['sequence_status'] = sequence_statuses[locus]
        row['model_sequence_length_bp'] = inference_metadata['sequence_length_bp']
        nearest = _nearest(ref_tss.get((str(row["gene_id"]), chrom, strand), []), summit)
        row["nearest_known_tss_distance_bp"] = nearest
        row["known_within_50bp"] = nearest is not None and nearest <= 50
        row["deeptssfilter_score"] = score
        row["filter_status"] = ('unscorable' if score is None else 'pass' if score >= threshold else 'fail')
        for position, molecules in sorted(row["_ctss"].items()):
            ctss_rows.append({"cluster_id": row["cluster_id"], "gene_id": row["gene_id"], "chrom": chrom,
                              "ctss_1": position, "strand": strand, "corrected_molecule_count": molecules})
    conn.close()
    _annotate_gene_overlap(raw_clusters, gene_intervals)
    high_confidence_clusters_path = (
        clusters_dir / f"{tso_mode}_route_final_high_confidence_clusters.tsv"
    )
    high_confidence_ctss_path = (
        ctss_dir / f"{tso_mode}_route_final_high_confidence_ctss.tsv"
    )
    final_raw_ctss_path = ctss_dir / f"{tso_mode}_route_final_raw_ctss.tsv"
    output_paths["final_raw_ctss"] = final_raw_ctss_path
    output_paths["final_high_confidence_clusters"] = high_confidence_clusters_path
    output_paths["final_high_confidence_ctss"] = high_confidence_ctss_path
    _write_tsv(output_paths["final_raw_clusters"], CLUSTER_COLUMNS, raw_clusters)
    _write_tsv(high_confidence_clusters_path, CLUSTER_COLUMNS,
               (row for row in raw_clusters if row["filter_status"] == "pass"))
    _write_tsv(final_raw_ctss_path, CTSS_COLUMNS, ctss_rows)
    passed = {row["cluster_id"] for row in raw_clusters if row["filter_status"] == "pass"}
    _write_tsv(high_confidence_ctss_path, CTSS_COLUMNS,
               (row for row in ctss_rows if row["cluster_id"] in passed))
    # Import lazily so reference/preflight commands do not pay AnnData's
    # import cost. Matrix reconstruction uses the retained evidence database,
    # so it exactly preserves cell barcodes and per-cell UMI correction.
    from .matrix import export_filtered_h5ad

    h5ad = export_filtered_h5ad(
        evidence_db=db_path,
        clusters_tsv=high_confidence_clusters_path,
        genes_tsv=Path(reference_tss_path).resolve().parent / "genes.tsv",
        outdir=root,
        umi_corrector=_directional_umi_correction,
        barcode_whitelist=barcode_whitelist,
        sample_id=sample_id,
        run_metadata={'model_sha256':inference_metadata.get('model_sha256',''),
                      'model_sequence_length_bp':inference_metadata['sequence_length_bp'],
                      'model_threshold':threshold,'tso_mode':tso_mode,
                      'endpoint_rule':resolved_endpoint_rule},
    )
    return {
        "cluster_splitting": asdict(split_config),
        "evidence_database": str(db_path), "raw_clusters": len(raw_clusters),
        "high_confidence_clusters": len(passed), "raw_ctss": len(ctss_rows),
        "output_files": {key: str(path) for key, path in output_paths.items()},
        "tso_mode": tso_mode,
        "endpoint_rule": resolved_endpoint_rule, "requested_endpoint_rule": endpoint_rule,
        "tso_sequence": tso_sequence if tso_mode in {"tso", "mix"} and tso_source == "evidence" else None,
        "assay_type": assay_type, "alignment_layout": resolved_layout,
        "evidence_mate": resolved_mate, "library_strand": library_strand,
        "tso_source": effective_tso_source,
        "configured_tso_source": tso_source,
        "tso_tag": tso_tag if effective_tso_source == "tag" else None,
        "evidence_technical_prefix_bp": evidence_technical_prefix_bp,
        "tso_junction_slop_bp": tso_junction_slop_bp,
        "tso_min_adjacent_match_bp": tso_min_adjacent_match_bp,
        "tso_reject_genomic_match": tso_reject_genomic_match,
        "raw_only": False, "min_gene_reads": min_gene_reads, "min_cluster_reads": min_cluster_reads,
        "off_known_tss_bp": off_known_tss_bp if tso_mode in {"off", "mix"} else None,
        "off_near_tss_bp": off_near_tss_bp if tso_mode in {"off", "mix"} else None,
        "off_proximal_tss_bp": off_proximal_tss_bp if tso_mode in {"off", "mix"} else None,
        "h5ad": h5ad,
        "model_inference": inference_metadata,
        "model_threshold": threshold,
        **dict(counters),
    }
