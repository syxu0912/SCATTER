"""Check that endpoint choice reaches stored evidence and final raw clusters."""
import csv
import sqlite3

import pysam
import pytest

from scatter.calling import call_clusters
from scatter.cli import make_parser


@pytest.mark.parametrize("rna_strand,flag,expected", [("+", 16, 1001), ("-", 0, 1150)])
def test_opposite_r2_boundary_is_used_before_clustering(tmp_path, rna_strand, flag, expected):
    reference = tmp_path / "reference_tss.tsv"
    reference.write_text("reference_tss_id\tgene_id\tchrom\ttss_1\tstrand\ttranscript_ids\n"
                         f"ref\tG\tchr1\t{expected}\t{rna_strand}\ttx\n")
    (tmp_path / "genes.tsv").write_text(
        "gene_id\tgene_name\tgene_type\tchrom\tstart_1\tend_1\tstrand\n"
        f"G\tG\tprotein_coding\tchr1\t1\t2000\t{rna_strand}\n")
    bam_path = tmp_path / "in.bam"
    with pysam.AlignmentFile(bam_path, "wb", header={"HD": {"SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 2000}]}) as bam:
        for i in range(60):
            read = pysam.AlignedSegment()
            read.query_name = f"m{i}"
            read.query_sequence = "A" * 150
            read.flag = flag
            read.reference_id = 0
            read.reference_start = 1000
            read.mapping_quality = 255
            read.cigartuples = [(0, 150)]
            read.set_tags([("CB", f"cell{i}"), ("UB", "AAAA"), ("GX", "G")])
            bam.write(read)
    pysam.index(str(bam_path))
    results = {}
    for rule in ("auto", "physical-5prime"):
        result = call_clusters(bam_path=str(bam_path), bam_index=str(bam_path)+".bai",
            fasta_path="unused.fa", reference_tss_path=str(reference), model_path="unused.pth",
            outdir=str(tmp_path / rule), alignment_layout="single", assay_type="5prime",
            evidence_mate="single", library_strand="gene", endpoint_rule=rule,
            tso_source="none", tso_mode="off", raw_only=True)
        assert result["eligible_alignments"] == 60
        with open(result["output_files"]["final_raw_clusters"]) as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        assert len(rows) == 1
        assert rows[0]["strand"] == rna_strand
        assert int(rows[0]["qualified_read_count"]) == 60
        with sqlite3.connect(result["evidence_database"]) as conn:
            assert {p for (p,) in conn.execute("SELECT DISTINCT pos FROM evidence")} == {int(rows[0]["summit_1"])}
        results[rule] = rows[0]
    assert int(results["auto"]["summit_1"]) == expected
    assert results["auto"]["endpoint_rule"] == "selected_alignment_transcript_5prime"
    assert int(results["physical-5prime"]["summit_1"]) == (1150 if rna_strand == "+" else 1001)
    assert results["physical-5prime"]["endpoint_rule"] == "selected_read_physical_5prime"


def test_cli_can_reproduce_old_endpoint_rule():
    args = make_parser().parse_args(["call", "--bam", "in.bam", "--gtf", "in.gtf",
        "--fasta", "in.fa", "--alignment-layout", "single", "--reference-tss", "ref.tsv",
        "--model", "model.pth", "--outdir", "out", "--endpoint-rule", "physical-5prime"])
    assert args.endpoint_rule == "physical-5prime"
