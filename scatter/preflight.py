"""Read-only validation for SCATTER inputs.

This module deliberately never changes a source BAM, FASTA, or GTF.  An index
may be created only when the caller explicitly supplies a destination path.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import pysam
from .io import open_text

from .protocol import (
    LIBRARY_STRANDS,
    PAIRED_R1_MIN_BP,
    TSO_SOURCES,
    alignment_strand,
    boolean_tag,
    clean_gene,
    forced_alignment_layout_for_r1,
    load_gene_catalog,
    record_matches_mate,
    resolve_alignment_layout,
    resolve_evidence_mate,
    transcript_strand,
)


SUPPORTED_LIBRARY_TYPES = {"10x_5prime_pe", "10x_5prime_se"}
TAG_PAIRS = (("CB", "UB"), ("CR", "UR"))


class PreflightError(RuntimeError):
    """Raised when an input violates the SCATTER input contract."""


@dataclass
class BamSummary:
    path: str
    index_path: str
    sort_order: str | None
    reference_count: int
    sampled_records: int
    paired_records: int
    proper_pair_records: int
    read1_records: int
    read2_records: int
    r1_query_length_histogram: dict[str, int]
    dominant_r1_query_length_bp: int | None
    mapped_records: int
    primary_records: int
    secondary_records: int
    supplementary_records: int
    mapq_histogram: dict[str, int]
    selected_cb_tag: str
    selected_umi_tag: str
    complete_tag_records: int
    tags_present: dict[str, int]
    gene_tag_records: int
    alignment_layout: str
    evidence_mate: str
    selected_mate_primary_mapped_records: int
    qualified_evidence_records: int
    alignment_same_as_gene_records: int
    alignment_opposite_gene_records: int
    resolved_strand_matches_gene_records: int
    library_strand: str
    tso_source: str
    tso_tag_records: int


@dataclass
class FastaSummary:
    path: str
    reference_count: int
    references: list[str]


@dataclass
class GtfSummary:
    path: str
    sampled_lines: int
    feature_counts: dict[str, int]
    has_gene_id: bool
    has_transcript_id: bool
    chromosomes_seen: list[str]


def _require_file(path: str | Path, label: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise PreflightError(f"{label} does not exist or is not a file: {resolved}")
    if not resolved.stat().st_size:
        raise PreflightError(f"{label} is empty: {resolved}")
    return resolved


def resolve_index_path(bam_path: str | Path, index_path: str | Path | None) -> Path:
    bam = _require_file(bam_path, "BAM")
    if index_path is not None:
        return _require_file(index_path, "BAM index")

    candidates = [
        Path(f"{bam}.bai"),
        bam.with_suffix(".bai"),
        Path(f"{bam}.csi"),
        bam.with_suffix(".csi"),
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    raise PreflightError(
        "No BAM index was found. Supply --bam-index, or create an index in the "
        "SCATTER output directory; source files are never modified by SCATTER."
    )


def select_tag_pair(
    tag_counts: Counter[str], cb_tag: str | None, umi_tag: str | None, whitelist: str | None
) -> tuple[str, str]:
    if (cb_tag is None) != (umi_tag is None):
        raise PreflightError("Specify --cb-tag and --umi-tag together.")
    if cb_tag is not None:
        pair = (cb_tag, umi_tag)
    elif tag_counts["CB"] and tag_counts["UB"]:
        pair = ("CB", "UB")
    elif tag_counts["CR"] and tag_counts["UR"]:
        pair = ("CR", "UR")
    else:
        raise PreflightError("No complete CB/UB or CR/UR tag pair was detected.")

    if pair == ("CR", "UR") and whitelist is None:
        raise PreflightError("CR/UR mode requires --barcode-whitelist.")
    return pair


def inspect_bam(
    bam_path: str | Path,
    index_path: str | Path,
    library_type: str | None,
    cb_tag: str | None = None,
    umi_tag: str | None = None,
    whitelist: str | None = None,
    sample_records: int = 100_000,
    *,
    alignment_layout: str | None = None,
    evidence_mate: str | None = "auto",
    reference_tss_path: str | Path | None = None,
    library_strand: str = "gene",
    tso_source: str = "evidence",
    tso_tag: str = "ZT",
    min_mapq: int = 30,
) -> BamSummary:
    try:
        resolved_layout = resolve_alignment_layout(library_type, alignment_layout)
        resolved_mate = resolve_evidence_mate(resolved_layout, evidence_mate)
    except ValueError as exc:
        raise PreflightError(str(exc)) from exc
    if library_strand not in LIBRARY_STRANDS:
        raise PreflightError(f"Unsupported --library-strand: {library_strand}")
    if tso_source not in TSO_SOURCES:
        raise PreflightError(f"Unsupported --tso-source: {tso_source}")
    if tso_source == "tag" and len(tso_tag) != 2:
        raise PreflightError("BAM TSO tag must be exactly two characters")
    gene_catalog = load_gene_catalog(reference_tss_path) if reference_tss_path is not None else None
    bam = _require_file(bam_path, "BAM")
    index = _require_file(index_path, "BAM index")
    if whitelist is not None:
        _require_file(whitelist, "Barcode whitelist")

    try:
        handle = pysam.AlignmentFile(str(bam), "rb", index_filename=str(index))
    except (OSError, ValueError) as exc:
        raise PreflightError(f"BAM/index cannot be opened together: {exc}") from exc

    with handle:
        header = handle.header.to_dict()
        sort_order = header.get("HD", {}).get("SO")
        if sort_order != "coordinate":
            raise PreflightError(f"BAM must be coordinate sorted; found SO={sort_order!r}.")

        counters: Counter[str] = Counter()
        tags: Counter[str] = Counter()
        mapqs: Counter[int] = Counter()
        complete_by_pair: Counter[tuple[str, str]] = Counter()
        r1_query_lengths: Counter[int] = Counter()
        pair_metrics: dict[tuple[str, str], Counter[str]] = {
            pair: Counter() for pair in TAG_PAIRS
        }
        if cb_tag is not None and umi_tag is not None:
            pair_metrics.setdefault((cb_tag, umi_tag), Counter())
        sampled = 0
        for record in handle.fetch(until_eof=True):
            if sampled >= sample_records:
                break
            sampled += 1
            counters["paired"] += int(record.is_paired)
            counters["proper_pair"] += int(record.is_proper_pair)
            counters["read1"] += int(record.is_read1)
            counters["read2"] += int(record.is_read2)
            counters["mapped"] += int(not record.is_unmapped)
            counters["secondary"] += int(record.is_secondary)
            counters["supplementary"] += int(record.is_supplementary)
            counters["primary"] += int(not record.is_secondary and not record.is_supplementary)
            mapqs[min(record.mapping_quality, 255)] += 1
            if record.is_read1 and not record.is_secondary and not record.is_supplementary:
                query_length = record.query_length
                if query_length:
                    # query_length includes soft clips, so STAR's
                    # --clip5pNbases does not hide the sequenced R1 length.
                    r1_query_lengths[int(query_length)] += 1

            # Cell Ranger BAMs routinely contain unmapped, secondary, or
            # supplementary records without the per-molecule CB/UB tags.
            # SCATTER discards those records before endpoint extraction, so
            # the tag contract must be checked on the same eligible class of
            # records rather than on every physical BAM record.
            if record.is_unmapped or record.is_secondary or record.is_supplementary:
                continue
            if not record_matches_mate(record, resolved_mate):
                continue
            counters["selected_mate_primary_mapped"] += 1
            counters["tag_check_records"] += 1
            record_tags = {tag for tag, _ in record.get_tags()}
            for tag in record_tags:
                tags[tag] += 1
            for pair in TAG_PAIRS:
                if set(pair).issubset(record_tags):
                    complete_by_pair[pair] += 1
            counters["gene_tag"] += int("GX" in record_tags or "GN" in record_tags)
            counters["tso_tag"] += int(tso_tag in record_tags)

            if record.mapping_quality < min_mapq or "GX" not in record_tags:
                continue
            gene = clean_gene(record.get_tag("GX"))
            if gene is None:
                continue
            gene_record = gene_catalog.get(gene) if gene_catalog is not None else None
            reference_name = handle.get_reference_name(record.reference_id)
            for pair, metrics in pair_metrics.items():
                if not set(pair).issubset(record_tags):
                    continue
                selected_cb_value = record.get_tag(pair[0])
                selected_umi_value = record.get_tag(pair[1])
                if (not selected_cb_value or selected_cb_value == "-" or
                        not selected_umi_value or selected_umi_value == "-"):
                    continue
                metrics["qualified_evidence"] += 1
                if tso_source == "tag":
                    try:
                        boolean_tag(record, tso_tag)
                        metrics["qualified_tso_tag"] += 1
                    except (KeyError, ValueError):
                        pass
                if gene_record is not None and gene_record.chrom == reference_name:
                    observed = alignment_strand(record)
                    resolved = transcript_strand(record, gene_record, library_strand)
                    metrics["gene_strand_known"] += 1
                    metrics["alignment_same_gene"] += int(observed == gene_record.strand)
                    metrics["alignment_opposite_gene"] += int(observed != gene_record.strand)
                    metrics["resolved_same_gene"] += int(resolved == gene_record.strand)

        if not sampled:
            raise PreflightError("BAM contains no records.")
        if resolved_layout == "paired" and counters["paired"] != sampled:
            raise PreflightError("Paired layout requires paired-end alignment flags.")
        if resolved_layout == "single" and counters["paired"]:
            raise PreflightError("Single layout requires single-end alignment flags.")
        if not counters["selected_mate_primary_mapped"]:
            raise PreflightError(
                f"Selected evidence mate {resolved_mate} has no primary mapped records in the preflight sample."
            )
        dominant_r1_length: int | None = None
        if resolved_layout == "paired" and resolved_mate == "read1":
            if not r1_query_lengths:
                raise PreflightError(
                    "Paired layout has no primary Read 1 sequence from which to enforce "
                    f"the {PAIRED_R1_MIN_BP}-bp R1 policy."
                )
            # Read lengths should normally be fixed. In mixed/trimmed input,
            # use the most frequent length and break exact ties conservatively
            # in favour of the shorter length.
            dominant_r1_length = max(
                r1_query_lengths, key=lambda length: (r1_query_lengths[length], -length)
            )
            if forced_alignment_layout_for_r1(dominant_r1_length) != "paired":
                raise PreflightError(
                    f"Dominant R1 query length is {dominant_r1_length} bp, below the forced "
                    f"paired threshold of {PAIRED_R1_MIN_BP} bp. Re-align this dataset in "
                    "unpaired/R2-only mode and use --alignment-layout single."
                )

        selected_cb, selected_umi = select_tag_pair(tags, cb_tag, umi_tag, whitelist)
        complete_selected = complete_by_pair[(selected_cb, selected_umi)]
        selected_metrics = pair_metrics[(selected_cb, selected_umi)]
        for metric in (
            "qualified_evidence", "qualified_tso_tag", "gene_strand_known",
            "alignment_same_gene", "alignment_opposite_gene", "resolved_same_gene",
        ):
            counters[metric] = selected_metrics[metric]
        tag_check_records = counters["tag_check_records"]
        if not complete_selected or not tag_check_records or complete_selected / tag_check_records < 0.5:
            raise PreflightError(
                f"Too few selected {resolved_mate} primary records contain {selected_cb}/{selected_umi}: "
                f"{complete_selected}/{tag_check_records} ({complete_selected / tag_check_records:.1%})."
            )
        if not counters["qualified_evidence"]:
            raise PreflightError(
                f"Selected evidence mate {resolved_mate} has no records passing MAPQ/GX/CB/UB checks "
                "in the preflight sample."
            )
        if tso_source == "tag" and counters["qualified_tso_tag"] != counters["qualified_evidence"]:
            raise PreflightError(
                f"Transferred TSO tag {tso_tag} is not valid and complete on qualified evidence records: "
                f"{counters['qualified_tso_tag']}/{counters['qualified_evidence']}."
            )
        if library_strand != "gene" and counters["gene_strand_known"] >= 100:
            agreement = counters["resolved_same_gene"] / counters["gene_strand_known"]
            if agreement < 0.95:
                raise PreflightError(
                    f"Configured library strand {library_strand} agrees with GX genes for only "
                    f"{agreement:.1%} of audited records; expected at least 95%."
                )
        return BamSummary(
            path=str(bam),
            index_path=str(index),
            sort_order=sort_order,
            reference_count=len(handle.references),
            sampled_records=sampled,
            paired_records=counters["paired"],
            proper_pair_records=counters["proper_pair"],
            read1_records=counters["read1"],
            read2_records=counters["read2"],
            r1_query_length_histogram={
                str(key): value for key, value in sorted(r1_query_lengths.items())
            },
            dominant_r1_query_length_bp=dominant_r1_length,
            mapped_records=counters["mapped"],
            primary_records=counters["primary"],
            secondary_records=counters["secondary"],
            supplementary_records=counters["supplementary"],
            mapq_histogram={str(key): value for key, value in sorted(mapqs.items())},
            selected_cb_tag=selected_cb,
            selected_umi_tag=selected_umi,
            complete_tag_records=complete_selected,
            tags_present=dict(sorted(tags.items())),
            gene_tag_records=counters["gene_tag"],
            alignment_layout=resolved_layout,
            evidence_mate=resolved_mate,
            selected_mate_primary_mapped_records=counters["selected_mate_primary_mapped"],
            qualified_evidence_records=counters["qualified_evidence"],
            alignment_same_as_gene_records=counters["alignment_same_gene"],
            alignment_opposite_gene_records=counters["alignment_opposite_gene"],
            resolved_strand_matches_gene_records=counters["resolved_same_gene"],
            library_strand=library_strand,
            tso_source=tso_source,
            tso_tag_records=counters["tso_tag"],
        )


def inspect_fasta(fasta_path: str | Path) -> FastaSummary:
    _require_file(fasta_path, "FASTA")
    # Preserve a run-local symlink's lexical path so its run-local .fai is
    # used. Resolving it would make htslib create an index beside the source.
    fasta = Path(fasta_path).expanduser().absolute()
    if not Path(str(fasta)+'.fai').is_file():
        raise PreflightError('FASTA requires a .fai index; scatter run can create one inside its output directory')
    try:
        with pysam.FastaFile(str(fasta)) as handle:
            references = list(handle.references)
    except (OSError, ValueError) as exc:
        raise PreflightError(
            f"FASTA requires a readable .fai index and could not be opened: {exc}"
        ) from exc
    return FastaSummary(path=str(fasta), reference_count=len(references), references=references)


def inspect_gtf(gtf_path: str | Path, sample_lines: int = 250_000) -> GtfSummary:
    gtf = _require_file(gtf_path, "GTF")
    features: Counter[str] = Counter()
    chromosomes: set[str] = set()
    has_gene_id = False
    has_transcript_id = False
    visited = 0
    with open_text(gtf) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            visited += 1
            chromosomes.add(fields[0])
            features[fields[2]] += 1
            attrs = fields[8]
            has_gene_id = has_gene_id or 'gene_id "' in attrs
            has_transcript_id = has_transcript_id or 'transcript_id "' in attrs
            if visited >= sample_lines:
                break
    required = {"gene", "transcript", "exon"}
    missing = sorted(required - set(features))
    if missing:
        raise PreflightError(f"GTF sample is missing required feature type(s): {', '.join(missing)}")
    if not has_gene_id or not has_transcript_id:
        raise PreflightError("GTF sample must contain gene_id and transcript_id attributes.")
    return GtfSummary(
        path=str(gtf),
        sampled_lines=visited,
        feature_counts=dict(sorted(features.items())),
        has_gene_id=has_gene_id,
        has_transcript_id=has_transcript_id,
        chromosomes_seen=sorted(chromosomes),
    )


def normalize_contig(name: str) -> str:
    token = name.lower()
    if token.startswith("chr"):
        token = token[3:]
    return {"m": "mt", "chrm": "mt"}.get(token, token)


def validate_contig_overlap(bam_path: str | Path, fasta: FastaSummary, gtf: GtfSummary) -> dict[str, list[str]]:
    bam = _require_file(bam_path, "BAM")
    with pysam.AlignmentFile(str(bam), "rb", check_sq=False) as handle:
        bam_refs = list(handle.references)
    fasta_normalized = {normalize_contig(item): item for item in fasta.references}
    gtf_normalized = {normalize_contig(item): item for item in gtf.chromosomes_seen}
    bam_normalized = {normalize_contig(item): item for item in bam_refs}
    missing_from_fasta = sorted(set(bam_normalized) - set(fasta_normalized))
    return {
        "bam_not_in_fasta": missing_from_fasta,
        # GTF inspection is deliberately bounded for a fast preflight, so this
        # is provenance only and must not be interpreted as a full GTF contig
        # compatibility test.
        "gtf_contigs_seen_in_sample": sorted(gtf_normalized.values()),
    }


def as_jsonable(value: object) -> object:
    if hasattr(value, "__dataclass_fields__"):
        return asdict(value)
    return value
