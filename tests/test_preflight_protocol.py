from __future__ import annotations

from pathlib import Path

import pysam
import pytest

from scatter.preflight import PreflightError, inspect_bam


def _write_reference(root: Path) -> Path:
    root.mkdir()
    reference_tss = root / "reference_tss.tsv"
    reference_tss.write_text("gene_id\tchrom\tstrand\ttss_1\nG1\tchr1\t+\t101\n")
    (root / "genes.tsv").write_text(
        "gene_id\tgene_name\tchrom\tstart_1\tend_1\tstrand\n"
        "G1\tG1\tchr1\t1\t10000\t+\n"
    )
    return reference_tss


def _segment(
    name: str, flag: int, start: int, umi: str, query_length: int = 100
) -> pysam.AlignedSegment:
    read = pysam.AlignedSegment()
    read.query_name = name
    read.flag = flag
    read.reference_id = 0
    read.reference_start = start
    read.mapping_quality = 60
    read.cigartuples = [(0, query_length)]
    read.query_sequence = "A" * query_length
    read.query_qualities = pysam.qualitystring_to_array("I" * query_length)
    read.set_tag("CB", "CELL-1")
    read.set_tag("UB", umi)
    read.set_tag("GX", "G1")
    return read


def _write_paired_bam(
    path: Path, *, include_read1: bool = True, read1_length: int = 100
) -> Path:
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 20000}]}
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for number in range(120):
            start = 100 + number * 20
            if include_read1:
                # Read 1 is opposite the + RNA strand: reverse-library rule.
                bam.write(_segment(
                    f"q{number}", 0x1 | 0x40 | 0x10, start,
                    f"A{number:03d}", read1_length,
                ))
            # Read 2 is on the + RNA strand: forward-library rule.
            bam.write(_segment(f"q{number}", 0x1 | 0x80, start + 10, f"B{number:03d}"))
    pysam.index(str(path))
    return Path(f"{path}.bai")


def test_preflight_accepts_explicit_mate_and_strand_conversion(tmp_path) -> None:
    reference_tss = _write_reference(tmp_path / "reference")
    bam = tmp_path / "paired.bam"
    bai = _write_paired_bam(bam)

    read1 = inspect_bam(
        bam, bai, None, sample_records=1_000, alignment_layout="paired",
        evidence_mate="read1", reference_tss_path=reference_tss,
        library_strand="reverse", tso_source="none",
    )
    read2 = inspect_bam(
        bam, bai, None, sample_records=1_000, alignment_layout="paired",
        evidence_mate="read2", reference_tss_path=reference_tss,
        library_strand="forward", tso_source="none",
    )

    assert read1.qualified_evidence_records == 120
    assert read1.dominant_r1_query_length_bp == 100
    assert read1.r1_query_length_histogram == {"100": 120}
    assert read1.resolved_strand_matches_gene_records == 120
    assert read2.qualified_evidence_records == 120
    assert read2.resolved_strand_matches_gene_records == 120


def test_preflight_forces_sub_100_bp_r1_to_unpaired(tmp_path) -> None:
    reference_tss = _write_reference(tmp_path / "reference")
    bam = tmp_path / "short_r1_paired.bam"
    bai = _write_paired_bam(bam, read1_length=99)

    with pytest.raises(PreflightError, match="Dominant R1 query length is 99 bp"):
        inspect_bam(
            bam, bai, None, sample_records=1_000, alignment_layout="paired",
            evidence_mate="read1", reference_tss_path=reference_tss,
            library_strand="reverse", tso_source="none",
        )


def test_preflight_rejects_wrong_strand_rule(tmp_path) -> None:
    reference_tss = _write_reference(tmp_path / "reference")
    bam = tmp_path / "paired.bam"
    bai = _write_paired_bam(bam)

    with pytest.raises(PreflightError, match="agrees with GX genes"):
        inspect_bam(
            bam, bai, None, sample_records=1_000, alignment_layout="paired",
            evidence_mate="read1", reference_tss_path=reference_tss,
            library_strand="forward", tso_source="none",
        )


def test_preflight_rejects_a_selected_mate_absent_from_bam(tmp_path) -> None:
    reference_tss = _write_reference(tmp_path / "reference")
    bam = tmp_path / "read2_only.bam"
    bai = _write_paired_bam(bam, include_read1=False)

    with pytest.raises(PreflightError, match="read1 has no primary mapped records"):
        inspect_bam(
            bam, bai, None, sample_records=1_000, alignment_layout="paired",
            evidence_mate="read1", reference_tss_path=reference_tss,
            library_strand="gene", tso_source="none",
        )
