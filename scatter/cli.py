"""Command line entry points for SCATTER."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Iterable

import pysam
import yaml

from . import __version__
from .clustering import SplitConfig
from .preflight import (
    PreflightError,
    as_jsonable,
    inspect_bam,
    inspect_fasta,
    inspect_gtf,
    resolve_index_path,
    validate_contig_overlap,
)
from .rampage import RampageParseError, parse_high_peak_bed
from .reference import build_reference
from .calling import (
    DEFAULT_MIN_CLUSTER_READS,
    DEFAULT_MIN_GENE_READS,
    DEFAULT_TSO_SEQUENCE,
    DEFAULT_UTR5_FLANK_BP,
    call_clusters,
)
from .protocol import ASSAY_TYPES, EVIDENCE_MATES, LIBRARY_STRANDS, TSO_SOURCES


def sha256(path: str | Path, chunk_size: int = 4 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def fingerprint(path: str | Path, mode: str) -> dict[str, object]:
    """Return reproducibility metadata without forcing a 72 GB reread by default."""
    target = Path(path).expanduser().resolve()
    info = target.stat()
    result: dict[str, object] = {
        "path": str(target),
        "size_bytes": info.st_size,
        "mtime_ns": info.st_mtime_ns,
        "mode": mode,
    }
    if mode == "full":
        result["sha256"] = sha256(target)
        return result
    if mode != "fast":
        raise ValueError(f"Unknown checksum mode: {mode}")
    chunk_size = min(4 * 1024 * 1024, info.st_size)
    digest = hashlib.sha256()
    with target.open("rb") as handle:
        digest.update(handle.read(chunk_size))
        if info.st_size > chunk_size:
            handle.seek(max(0, info.st_size - chunk_size))
            digest.update(handle.read(chunk_size))
    result["edge_sha256"] = digest.hexdigest()
    return result


def _output_dir(path: str | Path) -> Path:
    destination = Path(path).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    return destination


def build_index(bam_path: str, output_index: str) -> None:
    bam = Path(bam_path).expanduser().resolve()
    index = Path(output_index).expanduser().resolve()
    if not bam.is_file():
        raise PreflightError(f"BAM does not exist: {bam}")
    index.parent.mkdir(parents=True, exist_ok=True)
    # The positional output form is supported consistently by the deployed
    # htslib/samtools version and writes only the explicit destination index.
    pysam.index(str(bam), str(index))


def run_preflight(args: argparse.Namespace) -> int:
    outdir = _output_dir(args.outdir)
    metadata_dir = outdir / "metadata"
    metadata_dir.mkdir(exist_ok=True)
    index = resolve_index_path(args.bam, args.bam_index)
    bam = inspect_bam(
        args.bam,
        index,
        args.library_type,
        args.cb_tag,
        args.umi_tag,
        args.barcode_whitelist,
        args.sample_records,
        alignment_layout=args.alignment_layout,
        evidence_mate=args.evidence_mate,
        reference_tss_path=args.reference_tss,
        library_strand=args.library_strand,
        tso_source=args.tso_source,
        tso_tag=args.tso_tag,
        min_mapq=args.min_mapq,
    )
    fasta = inspect_fasta(args.fasta)
    gtf = inspect_gtf(args.gtf, args.gtf_sample_lines)
    contigs = validate_contig_overlap(args.bam, fasta, gtf)
    rampage = None
    if args.rampage_bed:
        peaks = parse_high_peak_bed(args.rampage_bed)
        rampage = {
            "path": str(Path(args.rampage_bed).expanduser().resolve()),
            "peak_count": len(peaks),
            "coordinate_convention": "BED summit is 0-based; SCATTER benchmark summit is 1-based",
            "chromosomes": sorted({peak.chrom for peak in peaks}),
        }
    report = {
        "scatter_version": __version__,
        "status": "PASS",
        "library_type": args.library_type,
        "alignment_layout": bam.alignment_layout,
        "assay_type": args.assay_type,
        "evidence_mate": bam.evidence_mate,
        "library_strand": args.library_strand,
        "tso_source": args.tso_source,
        "bam": as_jsonable(bam),
        "fasta": as_jsonable(fasta),
        "gtf": as_jsonable(gtf),
        "contig_validation": contigs,
        "rampage": rampage,
    }
    # argparse stores the dispatch function in the namespace; it is not a
    # user-facing parameter and cannot be serialized into a provenance file.
    config = {key: value for key, value in vars(args).items() if key != "handler"}
    config["resolved_bam_index"] = str(index)
    (metadata_dir / "preflight_report.json").write_text(json.dumps(report, indent=2) + "\n")
    (metadata_dir / "run_config.yaml").write_text(yaml.safe_dump(config, sort_keys=True))
    checksums = {
        "bam": fingerprint(args.bam, args.checksum_mode),
        "bam_index": fingerprint(index, args.checksum_mode),
        "fasta": fingerprint(args.fasta, args.checksum_mode),
        "gtf": fingerprint(args.gtf, args.checksum_mode),
    }
    if args.rampage_bed:
        checksums["rampage_bed"] = fingerprint(args.rampage_bed, args.checksum_mode)
    (metadata_dir / "input_checksums.json").write_text(json.dumps(checksums, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


def run_call(args: argparse.Namespace) -> int:
    from .model import bundled_model
    index = resolve_index_path(args.bam, args.bam_index)
    # Fail fast on a wrong mate, strand rule, or transferred-TSO contract
    # before creating a potentially very large evidence database.
    audit = inspect_bam(
        args.bam,
        index,
        args.library_type,
        args.cb_tag,
        args.umi_tag,
        args.barcode_whitelist,
        args.protocol_check_records,
        alignment_layout=args.alignment_layout,
        evidence_mate=args.evidence_mate,
        reference_tss_path=args.reference_tss,
        library_strand=args.library_strand,
        tso_source=args.tso_source if args.tso_mode in {"tso", "mix"} else "none",
        tso_tag=args.tso_tag,
        min_mapq=args.min_mapq,
    )
    result = call_clusters(
        bam_path=args.bam,
        bam_index=str(index),
        fasta_path=args.fasta,
        reference_tss_path=args.reference_tss,
        model_path=args.model or str(bundled_model()[0]),
        outdir=args.outdir,
        library_type=args.library_type,
        alignment_layout=args.alignment_layout,
        assay_type=args.assay_type,
        evidence_mate=args.evidence_mate,
        library_strand=args.library_strand,
        tso_source=args.tso_source,
        tso_tag=args.tso_tag,
        cb_tag=args.cb_tag or audit.selected_cb_tag,
        endpoint_rule=args.endpoint_rule,
        umi_tag=args.umi_tag or audit.selected_umi_tag,
        barcode_whitelist=args.barcode_whitelist,
        min_mapq=args.min_mapq,
        min_gene_reads=args.min_gene_reads,
        min_cluster_reads=args.min_cluster_reads,
        threshold=args.deeptssfilter_threshold,
        tso_mode=args.tso_mode,
        tso_sequence=args.tso_sequence,
        tso_min_match=args.tso_min_match,
        tso_max_mismatches=args.tso_max_mismatches,
        evidence_technical_prefix_bp=args.evidence_technical_prefix_bp,
        tso_junction_slop_bp=args.tso_junction_slop_bp,
        tso_min_adjacent_match_bp=args.tso_min_adjacent_match_bp,
        tso_reject_genomic_match=args.tso_reject_genomic_match,
        utr5_bed=args.utr5_bed,
        utr5_flank_bp=args.utr5_flank_bp,
        off_known_tss_bp=args.off_known_tss_bp,
        off_near_tss_bp=args.off_near_tss_bp,
        off_proximal_tss_bp=args.off_proximal_tss_bp,
        contigs=set(args.contigs.split(",")) if args.contigs else None,
        max_records=args.max_records,
        raw_only=args.raw_only,
        device=args.device,
        batch_size=args.batch_size,
        threads=args.threads,
        split_config=SplitConfig(
            max_width_bp=args.cluster_max_width_bp,
            valley_ratio=args.cluster_valley_ratio,
            smoothing_sigma_bp=args.cluster_smoothing_sigma_bp,
            min_peak_distance_bp=args.cluster_min_peak_distance_bp,
        ),
    )
    metadata_dir = _output_dir(args.outdir) / "metadata"
    metadata_dir.mkdir(exist_ok=True)
    config = {key: value for key, value in vars(args).items() if key != "handler"}
    config["resolved_bam_index"] = str(index)
    (metadata_dir / "call_config.yaml").write_text(yaml.safe_dump(config, sort_keys=True))
    (metadata_dir / "call_summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


def run_export_h5ad(args: argparse.Namespace) -> int:
    """Add AnnData matrices to a completed SCATTER call directory."""
    from .calling import _directional_umi_correction
    from .matrix import export_filtered_h5ad

    source = Path(args.scatter_outdir).expanduser().resolve()
    destination = Path(args.outdir).expanduser().resolve() if args.outdir else source
    summary_path = source / "metadata" / "call_summary.json"
    clusters_tsv = source / "clusters" / "high_confidence_clusters.tsv"
    if summary_path.is_file():
        summary = json.loads(summary_path.read_text())
        recorded = summary.get("output_files", {}).get("final_high_confidence_clusters")
        if recorded:
            clusters_tsv = Path(recorded).expanduser().resolve()
    result = export_filtered_h5ad(
        evidence_db=source / "work" / "ctss_evidence.sqlite",
        clusters_tsv=clusters_tsv,
        genes_tsv=Path(args.reference_tss).expanduser().resolve().parent / "genes.tsv",
        outdir=destination,
        umi_corrector=_directional_umi_correction,
        barcode_whitelist=args.barcode_whitelist,
    )
    print(json.dumps(result, indent=2))
    return 0


def add_common_inputs(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--bam", required=True)
    parser.add_argument("--bam-index")
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--fasta", required=True)
    layout = parser.add_mutually_exclusive_group(required=True)
    layout.add_argument("--alignment-layout", type=str.lower, choices=["paired", "single"],
                        help="Explicit BAM alignment layout")
    layout.add_argument("--library-type", choices=["10x_5prime_pe", "10x_5prime_se"],
                        help="Legacy alias for paired/single layout")
    parser.add_argument("--assay-type", type=str.lower, choices=sorted(ASSAY_TYPES), default="5prime")
    parser.add_argument("--evidence-mate", type=str.lower, choices=sorted(EVIDENCE_MATES), default="auto",
                        help="BAM mate that contributes one coordinate observation")
    parser.add_argument("--library-strand", type=str.lower, choices=sorted(LIBRARY_STRANDS), default="gene",
                        help=("gene: use GX gene strand; forward: alignment strand is RNA strand; "
                              "reverse: flip alignment strand"))
    parser.add_argument("--tso-source", type=str.lower, choices=sorted(TSO_SOURCES), default="evidence",
                        help="Detect TSO on evidence read, use transferred BAM tag, or disable")
    parser.add_argument("--tso-tag", default="ZT", help="Boolean BAM tag used with --tso-source tag")
    parser.add_argument("--cb-tag")
    parser.add_argument("--umi-tag")
    parser.add_argument("--barcode-whitelist")
    parser.add_argument("--outdir", required=True)


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="scatter", description="SCATTER TSS discovery")
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)

    from .pipeline import run, doctor, demo
    run_parser = sub.add_parser('run', help='Complete BAM-to-cell-TSS/CTSS workflow with the bundled 1000 bp Shift model')
    run_parser.add_argument('--bam', required=True)
    run_parser.add_argument('--bam-index')
    run_parser.add_argument('--fasta', required=True, help='FASTA used for the BAM alignment; bundled model is hg38')
    run_parser.add_argument('--gtf', required=True, help='Matching annotation GTF or GTF.gz')
    run_parser.add_argument('--outdir', required=True, help='New or empty output directory')
    run_parser.add_argument('--mode', required=True, choices=['tso','off','mix'])
    run_parser.add_argument('--layout', choices=['auto','paired','r2only'], default='auto')
    run_parser.add_argument('--evidence-mate', choices=['auto','read1','read2','single'], default='auto')
    cells = run_parser.add_mutually_exclusive_group(required=True)
    cells.add_argument('--barcodes', help='Called-cell barcode list, plain TSV or TSV.gz; order and zero-count cells are retained')
    cells.add_argument('--all-barcodes', action='store_true', help='Explicitly include all qualified BAM barcodes instead of a called-cell list')
    run_parser.add_argument('--sample-id')
    run_parser.add_argument('--cb-tag')
    run_parser.add_argument('--umi-tag')
    run_parser.add_argument('--tso-source', choices=['auto','evidence','tag'], default='auto')
    run_parser.add_argument('--tso-tag', default='ZT')
    run_parser.add_argument('--tso-sequence', default=DEFAULT_TSO_SEQUENCE)
    run_parser.add_argument('--technical-prefix-bp', type=int, default=0)
    run_parser.add_argument('--model', help='Optional replacement checkpoint; defaults to the packaged Shift hg38 1000 bp sigma=10 model')
    run_parser.add_argument('--threshold', type=float, default=.9)
    run_parser.add_argument('--device', default='auto', help='auto, cpu, cuda, or cuda:N')
    run_parser.add_argument('--threads', type=int, default=4)
    run_parser.add_argument('--batch-size', type=int, default=256)
    run_parser.add_argument('--min-mapq', type=int, default=30)
    run_parser.add_argument('--min-gene-reads', type=int, default=DEFAULT_MIN_GENE_READS)
    run_parser.add_argument('--min-cluster-reads', type=int, default=DEFAULT_MIN_CLUSTER_READS)
    run_parser.add_argument('--cluster-max-width-bp', type=int, default=200)
    run_parser.add_argument('--protocol-check-records', type=int, default=1000000)
    run_parser.add_argument('--contigs', help='Pilot only: restrict to comma-separated contigs')
    run_parser.add_argument('--max-records', type=int, help='Pilot only: cap scanned BAM records')
    run_parser.set_defaults(handler=run)
    check = sub.add_parser('doctor', help='Check installed dependencies and bundled model checksum')
    check.set_defaults(handler=doctor)
    example = sub.add_parser('demo', help='Run the complete workflow on bundled miniature paired/R2-only BAM data')
    example.add_argument('--outdir', required=True)
    example.add_argument('--layout', choices=['paired','r2only'], default='paired')
    example.add_argument('--mode', choices=['tso','off','mix'], default='mix')
    example.set_defaults(handler=demo)

    index = sub.add_parser("index", help="Create a BAM index at an explicit output path")
    index.add_argument("--bam", required=True)
    index.add_argument("--output-index", required=True)
    index.set_defaults(handler=lambda args: (build_index(args.bam, args.output_index) or 0))

    reference = sub.add_parser("build-reference", help="Build strand-aware gene and TSS reference catalogs")
    reference.add_argument("--gtf", required=True)
    reference.add_argument("--outdir", required=True)
    reference.add_argument("--merge-distance", type=int, default=50)
    reference.set_defaults(
        handler=lambda args: (
            print(json.dumps(build_reference(args.gtf, args.outdir, args.merge_distance), indent=2)) or 0
        )
    )

    preflight = sub.add_parser("preflight", help="Validate inputs without changing source files")
    add_common_inputs(preflight)
    preflight.add_argument("--reference-tss", required=True,
                           help="reference_tss.tsv whose sibling genes.tsv supplies GX strands")
    preflight.add_argument("--min-mapq", type=int, default=30)
    preflight.add_argument("--rampage-bed")
    preflight.add_argument("--sample-records", type=int, default=100_000)
    preflight.add_argument("--gtf-sample-lines", type=int, default=250_000)
    preflight.add_argument("--checksum-mode", choices=["fast", "full"], default="fast")
    preflight.set_defaults(handler=run_preflight)

    call = sub.add_parser("call", help="Discover raw and DeepTSSFilter-passing TSS clusters")
    add_common_inputs(call)
    call.add_argument("--reference-tss", required=True, help="reference_tss.tsv made by build-reference")
    call.add_argument(
        "--model",
        help="Optional checkpoint; default is the bundled 1000 bp Shift model",
    )
    call.add_argument("--min-mapq", type=int, default=30)
    call.add_argument("--endpoint-rule", choices=["auto", "physical-5prime", "transcript-5prime"], default="auto",
                      help="auto: RNA-facing 5-prime alignment boundary for single 5-prime input; physical read 5-prime otherwise")
    call.add_argument("--protocol-check-records", type=int, default=1_000_000,
                      help="Physical BAM records sampled for fail-fast mate/strand/TSO validation")
    call.add_argument("--min-gene-reads", type=int, default=DEFAULT_MIN_GENE_READS,
                      help=("Skip an evidence-route gene with fewer qualifying reads "
                            f"(default: {DEFAULT_MIN_GENE_READS})"))
    call.add_argument("--min-cluster-reads", "--min-raw-umis", dest="min_cluster_reads", type=int,
                      default=DEFAULT_MIN_CLUSTER_READS,
                      help=("Drop a route cluster with fewer qualifying reads "
                            f"(default: {DEFAULT_MIN_CLUSTER_READS})"))
    call.add_argument("--cluster-max-width-bp", type=int, default=200,
                      help="Maximum child-cluster span; 0 restores legacy unsplit boundaries (default: 200)")
    call.add_argument("--cluster-valley-ratio", type=float, default=0.25,
                      help="Split supported modes when the valley is at most this fraction of both peaks")
    call.add_argument("--cluster-smoothing-sigma-bp", type=float, default=5.0)
    call.add_argument("--cluster-min-peak-distance-bp", type=int, default=50)
    call.add_argument(
        "--deeptssfilter-threshold", type=float, default=0.9,
        help="DeepTSSFilter probability cutoff (default: 0.9)",
    )
    call.add_argument(
        "--tso-mode", choices=["off", "tso", "mix"], default="off",
        help=("off: all qualified reads then nearest-reference-TSS summit filtering; "
              "tso: unrestricted TSO-positive clusters; mix: independently cluster both "
              "routes, rescue off clusters overlapping TSO, then union overlaps"),
    )
    call.add_argument("--tso-sequence", default=DEFAULT_TSO_SEQUENCE)
    call.add_argument("--tso-min-match", type=int, default=8,
                      help="Minimum TSO suffix length; full sequence is used whenever available (default: 8)")
    call.add_argument("--tso-max-mismatches", type=int, default=2)
    call.add_argument(
        "--evidence-technical-prefix-bp", type=int, default=0,
        help="Sequenced 5′ technical prefix retained in the evidence soft clip (default: 0)",
    )
    call.add_argument(
        "--tso-junction-slop-bp", type=int, default=2,
        help="Maximum TSO offset from the soft-clip/alignment junction (default: 2)",
    )
    call.add_argument(
        "--tso-min-adjacent-match-bp", type=int, default=6,
        help="Minimum aligned anchor immediately beside the 5′ soft clip (default: 6)",
    )
    call.add_argument(
        "--allow-genomic-tso-match", dest="tso_reject_genomic_match", action="store_false",
        default=True, help="Do not reject TSO-like sequence encoded in the adjacent genome",
    )
    call.add_argument("--utr5-bed", help="Deprecated compatibility option; off filtering now uses --reference-tss")
    call.add_argument("--utr5-flank-bp", type=int, default=DEFAULT_UTR5_FLANK_BP,
                      help="Deprecated compatibility option; no longer used for off filtering")
    call.add_argument("--off-known-tss-bp", type=int, default=100,
                      help="Maximum summit distance for known-TSS class (default: 100)")
    call.add_argument("--off-near-tss-bp", type=int, default=500,
                      help="Maximum summit distance for near-TSS-shift class (default: 500)")
    call.add_argument("--off-proximal-tss-bp", type=int, default=1000,
                      help="Maximum summit distance for retained proximal-aTSS class (default: 1000)")
    call.add_argument("--contigs", help="Comma-separated pilot contigs; omitted scans all contigs")
    call.add_argument("--max-records", type=int, help="Read cap for a diagnostic pilot; omitted scans all records")
    call.add_argument("--raw-only", action="store_true",
                      help="Stop after writing route/final raw clusters and subdivision audit; no CTSS export")
    call.add_argument('--device', default='auto')
    call.add_argument('--threads', type=int, default=4)
    call.add_argument('--batch-size', type=int, default=256)
    call.set_defaults(handler=run_call)

    export = sub.add_parser(
        "export-h5ad",
        help="Export a completed filtered call as CamoTSS-compatible AnnData matrices",
    )
    export.add_argument("--scatter-outdir", required=True, help="Existing SCATTER call output directory")
    export.add_argument("--reference-tss", required=True, help="reference_tss.tsv used for the call")
    export.add_argument("--outdir", help="Matrix destination; defaults to --scatter-outdir")
    export.add_argument('--barcode-whitelist', help='Optional ordered list, including cells with zero retained counts')
    export.set_defaults(handler=run_export_h5ad)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = make_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.handler(args))
    except (PreflightError, RampageParseError, OSError, ValueError) as exc:
        print(f"SCATTER error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
