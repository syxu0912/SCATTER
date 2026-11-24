from __future__ import annotations

import pysam
import pytest

from scatter.protocol import (
    PAIRED_R1_MIN_BP,
    GeneRecord,
    biological_five_softclip,
    boolean_tag,
    evidence_five_position,
    evidence_five_softclip,
    forced_alignment_layout_for_r1,
    has_five_softclip_alignment_anchor,
    record_matches_mate,
    resolve_alignment_layout,
    resolve_evidence_mate,
    resolve_endpoint_rule,
    transcript_five_position,
    transcript_strand,
)


def _read(*, flag: int = 0, start: int = 100) -> pysam.AlignedSegment:
    read = pysam.AlignedSegment()
    read.query_name = "molecule"
    read.flag = flag
    read.reference_id = 0
    read.reference_start = start
    read.mapping_quality = 60
    read.cigartuples = [(4, 3), (0, 10), (4, 2)]
    read.query_sequence = "AAACCCCCCCCCCGG"
    return read


def test_paired_evidence_mate_can_select_read1_or_read2() -> None:
    read1 = _read(flag=0x1 | 0x40)
    read2 = _read(flag=0x1 | 0x80)

    assert record_matches_mate(read1, "read1")
    assert not record_matches_mate(read1, "read2")
    assert record_matches_mate(read2, "read2")
    assert not record_matches_mate(read2, "read1")


def test_layout_and_mate_resolution_reject_incompatible_combinations() -> None:
    assert resolve_alignment_layout(None, "paired") == "paired"
    assert resolve_alignment_layout("10x_5prime_se", None) == "single"
    assert resolve_evidence_mate("paired", "auto") == "read1"
    assert resolve_evidence_mate("paired", "read2") == "read2"
    assert resolve_evidence_mate("single", "auto") == "single"

    with pytest.raises(ValueError, match="Conflicting layout"):
        resolve_alignment_layout("10x_5prime_se", "paired")
    with pytest.raises(ValueError, match="Paired layout"):
        resolve_evidence_mate("paired", "single")
    with pytest.raises(ValueError, match="Single layout"):
        resolve_evidence_mate("single", "read1")


def test_r1_length_forces_layout_at_the_100_bp_boundary() -> None:
    assert PAIRED_R1_MIN_BP == 100
    assert forced_alignment_layout_for_r1(26) == "single"
    assert forced_alignment_layout_for_r1(99) == "single"
    assert forced_alignment_layout_for_r1(100) == "paired"
    assert forced_alignment_layout_for_r1(150) == "paired"

    with pytest.raises(ValueError, match="positive"):
        forced_alignment_layout_for_r1(0)


def test_forward_and_reverse_library_rules_transform_alignment_strand() -> None:
    reverse_alignment = _read(flag=0x10)
    plus_gene = GeneRecord("G1", "chr1", "+")

    assert transcript_strand(reverse_alignment, plus_gene, "forward") == "-"
    assert transcript_strand(reverse_alignment, plus_gene, "reverse") == "+"
    assert transcript_strand(reverse_alignment, plus_gene, "gene") == "+"


def test_alignment_orientation_controls_physical_five_prime_evidence_end() -> None:
    forward = _read(start=100)
    reverse = _read(flag=0x10, start=100)

    assert evidence_five_position(forward) == 101
    assert evidence_five_position(reverse) == 110
    assert evidence_five_softclip(forward) == "AAA"
    assert evidence_five_softclip(reverse) == "GG"


def test_r2_opposite_strand_uses_rna_facing_boundary() -> None:
    reverse = _read(flag=0x10, start=1000)
    reverse.cigartuples = [(0, 150)]
    assert evidence_five_position(reverse) == 1150
    assert transcript_five_position(reverse, "+") == 1001
    forward = _read(start=1000)
    forward.cigartuples = [(0, 150)]
    assert evidence_five_position(forward) == 1001
    assert transcript_five_position(forward, "-") == 1150


def test_rna_boundary_respects_splicing_and_does_not_extrapolate_softclips() -> None:
    read = _read(flag=0x10, start=1000)
    read.cigartuples = [(4, 5), (0, 50), (3, 900), (0, 100), (4, 7)]
    assert transcript_five_position(read, "+") == 1001
    assert transcript_five_position(read, "-") == 2050
    with pytest.raises(ValueError, match="Invalid transcript strand"):
        transcript_five_position(read, ".")


def test_endpoint_default_and_explicit_legacy_rule() -> None:
    assert resolve_endpoint_rule("auto", "5prime", "single") == "transcript-5prime"
    for assay, layout in [("5prime", "paired"), ("3prime", "single"), ("full_length", "single")]:
        assert resolve_endpoint_rule("auto", assay, layout) == "physical-5prime"
    assert resolve_endpoint_rule("physical-5prime", "5prime", "single") == "physical-5prime"
    with pytest.raises(ValueError, match="Unknown endpoint rule"):
        resolve_endpoint_rule("invalid", "5prime", "single")


def test_technical_prefix_is_removed_from_the_orientation_aware_side() -> None:
    forward = _read(start=100)
    reverse = _read(flag=0x10, start=100)

    assert biological_five_softclip(forward, 2) == "A"
    assert biological_five_softclip(reverse, 1) == "G"
    assert biological_five_softclip(forward, 3) == ""
    assert has_five_softclip_alignment_anchor(forward, 6)
    assert has_five_softclip_alignment_anchor(reverse, 6)


@pytest.mark.parametrize(
    ("value", "expected"),
    [(1, True), (0, False), ("positive", True), ("negative", False)],
)
def test_boolean_tso_tag(value: object, expected: bool) -> None:
    read = _read()
    read.set_tag("ZT", value)
    assert boolean_tag(read, "ZT") is expected
