from __future__ import annotations

import csv
import json
import random
import sqlite3

import numpy as np
import pysam
import pytest
import yaml

from scatter.annotation import GeneIntervalIndex, annotate_gene_overlap
from scatter.calling import _route_candidates, call_clusters
from scatter.clustering import SplitConfig, split_sites
from scatter.cli import make_parser


def test_long_continuous_component_is_bounded_and_conserves_sites():
    sites = [(p, 5) for p in range(100, 2101, 10)]
    children = split_sites(sites, SplitConfig(), 50)
    assert [site for piece, _ in children for site in piece] == sites
    assert all(piece[-1][0] - piece[0][0] + 1 <= 200 for piece, _ in children)
    assert all("max_width" in reason for _, reason in children)


def test_supported_peaks_are_split_at_valley_even_below_width_cap():
    sites = [(p, 60 if p in (100, 200) else 1) for p in range(100, 201, 10)]
    children = split_sites(sites, SplitConfig(), 50)
    assert len(children) == 2
    assert all(reason == "valley" for _, reason in children)
    assert all(sum(c for _, c in piece) >= 50 for piece, _ in children)
    assert children[0][0][-1][0] < children[1][0][0][0]


def test_single_peak_and_exact_width_boundary_are_not_arbitrarily_split():
    sites = [(100, 1), (199, 100), (299, 1)]
    assert split_sites(sites, SplitConfig(), 50) == [(sites, "unchanged")]
    assert split_sites([(100, 100)], SplitConfig(max_width_bp=1), 50) == [([(100, 100)], "unchanged")]


def test_zero_width_restores_legacy_and_invalid_configuration_fails():
    sites = [(0, 100), (5000, 100)]
    assert split_sites(sites, SplitConfig(max_width_bp=0), 50) == [(sites, "legacy")]
    for kwargs in ({"max_width_bp": -1}, {"valley_ratio": 1}, {"smoothing_sigma_bp": 0}):
        with pytest.raises(ValueError):
            SplitConfig(**kwargs)


def test_split_site_conservation_for_varied_profiles():
    rng = random.Random(734)
    for _ in range(35):
        positions = sorted(rng.sample(range(1, 2500), rng.randrange(2, 90)))
        sites = [(p, rng.randrange(1, 100)) for p in positions]
        children = split_sites(sites, SplitConfig(max_width_bp=150), 50)
        assert [s for piece, _ in children for s in piece] == sites
        assert all(piece[-1][0] - piece[0][0] + 1 <= 150 for piece, _ in children)


def test_route_split_audit_accounts_for_children_below_original_read_threshold():
    conn = sqlite3.connect(":memory:")
    conn.execute("CREATE TABLE evidence(chrom,gene,strand,pos,cb,umi,n_reads,tso_reads)")
    # One legacy component held together by intermediate sites; most children
    # are too weak to pass 50, and must be counted in the discard audit.
    conn.executemany("INSERT INTO evidence VALUES(?,?,?,?,?,?,?,?)", [
        ("chr1", "G", "+", p, f"cell{p}", "AAAA", 80 if p == 100 else 1, 0)
        for p in range(100, 701, 30)
    ])
    audit = []
    rows = _route_candidates(conn, "off", split_config=SplitConfig(max_width_bp=100), split_audit=audit)
    assert len(audit) == 1
    assert audit[0]["discarded_reads"] > 0
    assert audit[0]["parent_reads"] == audit[0]["retained_reads"] + audit[0]["discarded_reads"]
    assert audit[0]["retained_reads"] == sum(r["qualified_read_count"] for r in rows)
    assert all(r["qualified_read_count"] >= 50 and r["end_1"]-r["start_1"]+1 <= 100 for r in rows)
    conn.close()


def test_interval_queries_do_not_inherit_genes_from_previous_wider_query():
    index = GeneIntervalIndex([(1, 1500, "A", "+"), (500, 600, "B", "-")])
    assert {x[2] for x in index.overlaps(10, 1000)} == {"A", "B"}
    assert {x[2] for x in index.overlaps(20, 30)} == {"A"}
    assert {x[2] for x in index.overlaps(600, 600)} == {"A", "B"}
    assert {x[2] for x in index.overlaps(601, 601)} == {"A"}


def test_interval_index_matches_bruteforce_for_nested_and_unordered_queries():
    rng = random.Random(907)
    intervals = []
    for i in range(250):
        a = rng.randrange(1, 2000)
        intervals.append((a, a+rng.randrange(1, 800), str(i), "+" if i % 2 else "-"))
    index = GeneIntervalIndex(list(intervals))
    for _ in range(180):
        a = rng.randrange(1, 2300)
        b = a+rng.randrange(0, 900)
        assert set(index.overlaps(a, b)) == {x for x in intervals if x[0] <= b and x[1] >= a}


def test_unique_gx_survives_antisense_and_nested_span_overlap():
    genes = {"chr1": GeneIntervalIndex([
        (1, 1000, "A", "+"), (100, 500, "antisense", "-"), (100, 300, "nested", "+"),
    ])}
    rows = [{"chrom": "chr1", "gene_id": "A", "strand": "+", "start_1": 150,
             "end_1": 200, "filter_status": "pass"}]
    annotate_gene_overlap(rows, genes)
    row = rows[0]
    assert row["filter_status"] == "pass" and row["gene_assignment_status"] == "unique"
    assert row["same_strand_overlapping_gene_ids"] == "A;nested"
    assert row["opposite_strand_overlapping_gene_ids"] == "antisense"


@pytest.mark.parametrize("chrom,gene,strand", [("chr2", "A", "+"), ("chr1", "unknown", "+"), ("chr1", "A", "-")])
def test_invalid_gx_reference_is_still_excluded(chrom, gene, strand):
    rows = [{"chrom": chrom, "gene_id": gene, "strand": strand, "start_1": 10,
             "end_1": 20, "filter_status": "pass"}]
    annotate_gene_overlap(rows, {"chr1": GeneIntervalIndex([(1, 100, "A", "+")])})
    assert rows[0]["filter_status"] == "fail_ambiguous_gene"


def test_gx_supported_endpoint_extension_is_flagged_without_reassignment():
    row = {"chrom": "chr1", "gene_id": "A", "strand": "+", "start_1": 50,
           "end_1": 70, "filter_status": "pass"}
    annotate_gene_overlap([row], {"chr1": GeneIntervalIndex([(100, 500, "A", "+")])})
    assert row["filter_status"] == "pass"
    assert row["assigned_gene_span_overlap"] is False
    assert row["gene_assignment_evidence"] == "unique_GX"


def test_cli_records_split_configuration():
    args = make_parser().parse_args([
        "call", "--bam", "in.bam", "--gtf", "in.gtf", "--fasta", "in.fa",
        "--alignment-layout", "single", "--reference-tss", "ref.tsv", "--model", "model.pth",
        "--outdir", "out", "--cluster-max-width-bp", "300",
    ])
    assert args.cluster_max_width_bp == 300
    assert args.cluster_valley_ratio == 0.25


def test_full_call_preserves_antisense_supported_gx_and_exports_conserved_matrix(tmp_path, monkeypatch):
    from scatter import calling
    import anndata as ad

    reference = tmp_path / "reference_tss.tsv"
    reference.write_text("reference_tss_id\tgene_id\tchrom\ttss_1\tstrand\ttranscript_ids\nref\tA\tchr1\t400\t+\ttx\n")
    (tmp_path / "genes.tsv").write_text(
        "gene_id\tgene_name\tgene_type\tchrom\tstart_1\tend_1\tstrand\n"
        "A\tA\tprotein_coding\tchr1\t1\t1200\t+\nB\tB\tlncRNA\tchr1\t1\t1200\t-\n"
    )
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\n" + "ACGT" * 400 + "\n")
    pysam.faidx(str(fasta))
    bam = tmp_path / "evidence.bam"
    with pysam.AlignmentFile(bam, "wb", header={"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1600}]}) as f:
        for i in range(120):
            read = pysam.AlignedSegment()
            read.query_name = f"r{i}"
            read.query_sequence = "A" * 30
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 399 if i < 60 else 459
            read.mapping_quality = 255
            read.cigar = [(0, 30)]
            read.set_tags([("CB", f"cell{i}"), ("UB", "AAAA"), ("GX", "A")])
            f.write(read)
    pysam.index(str(bam))
    from scatter import model as inference
    def scored(rows, *args, **kwargs):
        keys = [(r['chrom'], int(r['summit_1']), r['strand']) for r in rows]
        return ({k:.95 for k in keys}, {k:'valid' for k in keys}, {'sequence_length_bp':1000})
    monkeypatch.setattr(inference, 'score_loci', scored)
    result = call_clusters(bam_path=str(bam), bam_index=str(bam)+".bai", fasta_path=str(fasta),
                           reference_tss_path=str(reference), model_path="unused.pth", outdir=str(tmp_path/"out"),
                           alignment_layout="single", evidence_mate="single", tso_mode="off")
    assert result["high_confidence_clusters"] == 2
    matrix = ad.read_h5ad(result["h5ad"]["files"]["scTSS_count_all"])
    assert int(matrix.X.sum()) == 120
    assert set(matrix.var["opposite_strand_overlapping_gene_ids"]) == {"B"}
    with open(result["output_files"]["cluster_splitting_audit"]) as f:
        audit = list(csv.DictReader(f, delimiter="\t"))
    assert sum(int(r["retained_reads"]) + int(r["discarded_reads"]) for r in audit) == 120

    # Legacy parent catalogs + immutable evidence must produce the exact
    # same revised raw tables as a fresh call from the original small BAM.
    from scatter.recluster import recluster_run
    legacy_dir = tmp_path / "legacy"
    legacy = call_clusters(bam_path=str(bam), bam_index=str(bam)+".bai", fasta_path=str(fasta),
        reference_tss_path=str(reference), model_path="unused", outdir=str(legacy_dir),
        alignment_layout="single", tso_mode="off", raw_only=True,
        split_config=SplitConfig(max_width_bp=0))
    assert legacy["raw_clusters"] == 1
    (legacy_dir/"metadata/call_summary.json").write_text(json.dumps(legacy))
    (legacy_dir/"metadata/call_config.yaml").write_text(yaml.safe_dump({"reference_tss": str(reference)}))
    revised = recluster_run(legacy_dir, tmp_path/"revised")
    fresh = call_clusters(bam_path=str(bam), bam_index=str(bam)+".bai", fasta_path=str(fasta),
        reference_tss_path=str(reference), model_path="unused", outdir=str(tmp_path/"fresh"),
        alignment_layout="single", tso_mode="off", raw_only=True)
    with open(revised["output_files"]["final_raw_clusters"]) as f:
        reclustered_rows = list(csv.DictReader(f, delimiter="\t"))
    with open(fresh["output_files"]["final_raw_clusters"]) as f:
        fresh_rows = list(csv.DictReader(f, delimiter="\t"))
    assert reclustered_rows == fresh_rows
    assert revised["split_read_accounting"]["off"]["retained_reads"] == 120
    with pytest.raises(FileExistsError):
        recluster_run(legacy_dir, tmp_path/"revised")
