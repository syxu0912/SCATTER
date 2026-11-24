"""Recompute raw candidates from immutable evidence and legacy parent catalogs.

Legacy catalogs already contain every adaptive-gap component meeting the
original threshold. Subdivision cannot promote a smaller discarded component,
so they are sufficient to rerun the two fixes without rereading source BAMs.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import asdict
from datetime import datetime, timezone
import csv
import hashlib
import json
from pathlib import Path
import sqlite3
import time

import yaml

from . import __version__
from .annotation import load_genes, annotate_gene_overlap
from .clustering import SplitConfig, split_sites
from .calling import (
    CLUSTER_COLUMNS, SPLIT_AUDIT_COLUMNS, _add_protocol_metadata,
    _annotate_reference_tss_support, _assign_route_cluster_ids,
    _filter_off_candidates, _load_reference_tss, _merge_mix_candidates,
    _output_paths, _quantify_route_piece, _write_tsv,
    _route_piece_row, _weighted_median,
)


INTEGER_FIELDS = {
    "start_1", "end_1", "summit_1", "adaptive_gap_bp", "qualified_read_count",
    "raw_site_umi_count", "corrected_molecule_count", "ctss_count",
    "nearest_known_tss_distance_bp", "tso_matched_alignment_count", "tso_total_alignment_count",
    "tso_read_count", "off_read_count", "tso_cluster_count", "off_cluster_count",
    "source_cluster_count", "nearest_reference_tss_1", "signed_reference_tss_distance_bp",
    "utr5_flank_bp",
}


def _read_parents(path: Path) -> list[dict[str, object]]:
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        for key, value in list(row.items()):
            if key in INTEGER_FIELDS and value != "":
                row[key] = int(value)
            elif key == "tso_fraction_of_alignments" and value != "":
                row[key] = float(value)
    return rows


def _hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _status(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(payload | {"updated_utc": datetime.now(timezone.utc).isoformat()}, indent=2) + "\n")
    temporary.replace(path)


def _recluster_route(conn, parents, route, config, minimum, audit, progress):
    result = []
    expression = "tso_reads" if route == "tso" else "n_reads"
    total_expression = "CASE WHEN tso_reads>0 THEN n_reads ELSE 0 END" if route == "tso" else "n_reads"
    for number, parent in enumerate(parents, 1):
        chrom, gene, strand = (str(parent[k]) for k in ("chrom", "gene_id", "strand"))
        start, end = int(parent["start_1"]), int(parent["end_1"])
        site_totals = conn.execute(
            f"SELECT pos,SUM({expression}),SUM({total_expression}),SUM(tso_reads) FROM evidence "
            "WHERE chrom=? AND gene=? AND strand=? AND pos BETWEEN ? AND ? "
            f"GROUP BY pos HAVING SUM({expression})>0 ORDER BY pos",
            (chrom, gene, strand, start, end),
        ).fetchall()
        sites = [(p, selected) for p, selected, _, _ in site_totals]
        totals = {p: (n, tso) for p, _, n, tso in site_totals}
        parent_reads = sum(c for _, c in sites)
        if parent_reads != int(parent["qualified_read_count"]):
            raise ValueError(f"Source evidence disagrees with {route} parent {parent['cluster_id']}")
        unique_molecules = parent_reads == int(parent["corrected_molecule_count"])
        if unique_molecules:
            # Each corrected molecule contributes one while each record has
            # a positive integer weight. Equality proves every selected UMI
            # group has exactly one read; all child subsets inherit this.
            if (len(sites) != int(parent["ctss_count"])
                    or _weighted_median(dict(sites)) != int(parent["summit_1"])):
                raise ValueError("Source parent contradicts its unique-molecule count contract")
        children = split_sites(sites, config, minimum)
        kept, discarded, dropped = 0, 0, 0
        parent_id = f"{route.upper()}_PARENT_{number:08d}"
        for piece, reason in children:
            reads = sum(c for _, c in piece)
            if reads < minimum:
                discarded += reads
                dropped += 1
                continue
            if len(children) == 1:
                # Same interval and unchanged UMI/CTSS policy: existing exact
                # quantification is reusable. Source counts were checked above.
                row = dict(parent)
                for key in ("filter_status", "off_route_filter_status", "off_route_filter_reason"):
                    row.pop(key, None)
            else:
                if unique_molecules:
                    row = _route_piece_row(chrom, gene, strand, piece,
                        int(parent["adaptive_gap_bp"]), route, Counter(dict(piece)), reads,
                        sum(totals[p][0] for p, _ in piece), sum(totals[p][1] for p, _ in piece))
                else:
                    row = _quantify_route_piece(conn, chrom, gene, strand, piece,
                                                int(parent["adaptive_gap_bp"]), route)
                if row is None:
                    raise ValueError("A supported child has no quantified evidence")
                row.pop("_ctss", None)
            row.update(parent_cluster_id=parent_id, parent_start_1=start, parent_end_1=end,
                       parent_qualified_read_count=parent_reads, cluster_split_reason=reason,
                       cluster_max_width_bp=config.max_width_bp)
            result.append(row)
            kept += reads
        if kept + discarded != parent_reads:
            raise RuntimeError("Read conservation failed during reclustering")
        audit.append(dict(
            parent_cluster_id=parent_id, route=route, chrom=chrom, gene_id=gene, strand=strand,
            start_1=start, end_1=end, parent_width_bp=end-start+1, parent_reads=parent_reads,
            retained_children=sum(sum(c for _, c in p) >= minimum for p, _ in children),
            discarded_children=dropped, retained_reads=kept, discarded_reads=discarded,
            split_reason="+".join(sorted({reason for _, reason in children})),
            max_child_width_bp=max(p[-1][0]-p[0][0]+1 for p, _ in children),
        ))
        if number % 500 == 0 or number == len(parents):
            progress(route, number, len(parents), len(result))
    return result


def recluster_run(source_dir: str | Path, outdir: str | Path,
                  split_config: SplitConfig | None = None) -> dict:
    source, destination = Path(source_dir).resolve(), Path(outdir).resolve()
    config = split_config or SplitConfig()
    if destination.exists() and any(destination.iterdir()):
        raise FileExistsError(f"Refusing to reuse a nonempty output directory: {destination}")
    if source == destination or source in destination.parents:
        raise ValueError("Revised outputs must be outside the source run directory")
    summary_path, config_path = source/"metadata/call_summary.json", source/"metadata/call_config.yaml"
    previous = json.loads(summary_path.read_text())
    original_config = yaml.safe_load(config_path.read_text())
    if not previous.get("raw_only"):
        raise ValueError("Expected an immutable raw-only source run")
    mode = previous["tso_mode"]
    if mode not in {"off", "tso", "mix"}:
        raise ValueError(f"Unsupported route mode: {mode}")
    minimum = int(previous["min_cluster_reads"])
    evidence = Path(previous["evidence_database"]).resolve()
    before = evidence.stat()
    reference = Path(original_config["reference_tss"]).resolve()
    for folder in ("clusters", "metadata", "work"):
        (destination/folder).mkdir(parents=True, exist_ok=True)
    (destination/"work/ctss_evidence.sqlite").symlink_to(evidence)
    started = time.monotonic()
    status_path = destination/"metadata/progress.json"

    def progress(route, done, total, children):
        payload = dict(status="running", route=route, parents_done=done, parents_total=total,
                       retained_children=children, elapsed_seconds=round(time.monotonic()-started, 2))
        _status(status_path, payload)
        print(json.dumps(payload), flush=True)

    audit = []
    provenance = {}
    candidates = {"tso": [], "off": []}
    with sqlite3.connect(evidence.as_uri()+"?mode=ro&immutable=1", uri=True) as conn:
        for route, key in (("tso", "tso_route_raw_clusters"),
                           ("off", "off_route_before_reference_tss_filter_raw_clusters")):
            if key not in previous["output_files"]:
                continue
            path = Path(previous["output_files"][key]).resolve()
            parents = _read_parents(path)
            provenance[route] = {"path": str(path), "sha256": _hash(path), "parents": len(parents)}
            progress(route, 0, len(parents), 0)
            candidates[route] = _recluster_route(conn, parents, route, config, minimum, audit, progress)
    tso, off = candidates["tso"], candidates["off"]
    _annotate_reference_tss_support(tso+off, _load_reference_tss(reference),
        known_bp=int(previous["off_known_tss_bp"] or 100),
        near_bp=int(previous["off_near_tss_bp"] or 500),
        proximal_bp=int(previous["off_proximal_tss_bp"] or 1000))
    _assign_route_cluster_ids(tso, "TSO_RAW")
    _assign_route_cluster_ids(off, "OFF_PRETSS")
    reference_off = _filter_off_candidates([dict(row) for row in off], [])
    mix_off = _filter_off_candidates(off, tso) if mode == "mix" else reference_off
    final = tso if mode == "tso" else reference_off if mode == "off" else _merge_mix_candidates(tso, mix_off)
    all_rows = list({id(row): row for row in tso+off+reference_off+mix_off+final}.values())
    _add_protocol_metadata(all_rows, assay_type=previous["assay_type"],
        alignment_layout=previous["alignment_layout"], evidence_mate=previous["evidence_mate"],
        library_strand=previous["library_strand"], tso_source=previous["tso_source"],
        tso_tag=previous.get("tso_tag") or "ZT",
        endpoint_rule=previous.get("endpoint_rule", "physical-5prime"))
    annotate_gene_overlap(all_rows, load_genes(reference))
    paths = _output_paths(destination, mode)
    for row in off:
        row["tso_mode"] = mode
        row["filter_status"] = "pending_reference_tss_filter"
    for key, rows in (
        ("tso_route_raw_clusters", tso),
        ("off_route_before_reference_tss_filter_raw_clusters", off),
        ("off_route_after_reference_tss_filter_raw_clusters", reference_off),
        ("mix_route_off_after_tso_rescue_raw_clusters", mix_off),
    ):
        if key in paths:
            _write_tsv(paths[key], CLUSTER_COLUMNS, rows)
    for number, row in enumerate(final, 1):
        row.update(cluster_id=f"SCATTER_{number:08d}", tso_mode=mode, filter_status="raw_only")
    _write_tsv(paths["final_raw_clusters"], CLUSTER_COLUMNS, final)
    audit_path = destination/"metadata/cluster_splitting_audit.tsv"
    _write_tsv(audit_path, SPLIT_AUDIT_COLUMNS, audit)
    paths["cluster_splitting_audit"] = audit_path
    result = previous | dict(
        scatter_version=__version__, raw_clusters=len(final), raw_only=True,
        evidence_database=str(destination/"work/ctss_evidence.sqlite"),
        source_run=str(source), source_evidence_database=str(evidence),
        source_call_summary_sha256=_hash(summary_path), source_parent_catalogs=provenance,
        cluster_splitting=asdict(config), output_files={k: str(p) for k,p in paths.items()},
        off_route_raw_clusters_before_tss_filter=len(off),
        off_route_clusters_retained_by_reference_tss=len(reference_off),
        off_route_clusters_filtered_by_tss_distance=len(off)-len(reference_off),
        elapsed_seconds=round(time.monotonic()-started, 2),
    )
    if mode in {"tso", "mix"}:
        result["tso_route_raw_clusters"] = len(tso)
    if mode == "mix":
        result["mix_off_clusters_after_tso_rescue"] = len(mix_off)
        result["mix_off_clusters_rescued_by_tso_overlap"] = sum(r["off_route_filter_reason"] == "tso_overlap_rescue" for r in mix_off)
    split_stats = {}
    for route in candidates:
        records = [r for r in audit if r["route"] == route]
        if records:
            split_stats[route] = {k: sum(r[k] for r in records) for k in (
                "parent_reads", "retained_reads", "discarded_reads", "retained_children", "discarded_children")}
    result["split_read_accounting"] = split_stats
    result["final_source_counts"] = dict(Counter(r["candidate_source"] for r in final))
    widths = sorted(int(r["end_1"])-int(r["start_1"])+1 for r in final)
    result["final_widths"] = {"maximum_bp": max(widths, default=0),
        "median_bp": widths[len(widths)//2] if widths else 0,
        "over_1000bp": sum(w > 1000 for w in widths)}
    if config.max_width_bp and max(widths, default=0) > config.max_width_bp:
        raise RuntimeError("Final candidates exceeded the configured width cap")
    after = evidence.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise RuntimeError("Source evidence changed during read-only reprocessing")
    (destination/"metadata/call_summary.json").write_text(json.dumps(result, indent=2)+"\n")
    revised_config = original_config | dict(command="recluster-raw", outdir=str(destination),
        source_run=str(source), scatter_version=__version__, cluster_splitting=asdict(config))
    (destination/"metadata/call_config.yaml").write_text(yaml.safe_dump(revised_config, sort_keys=True))
    _status(status_path, {"status": "complete", "raw_clusters": len(final),
                          "elapsed_seconds": result["elapsed_seconds"]})
    return result
