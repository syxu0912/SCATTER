from __future__ import annotations

import pysam

from scatter.calling import _has_tso, _reference_matches_tso, _tso_match_length
from scatter.deeptssfilter import reverse_complement


TSO = "TTTCTTATATGGG"
PREFIX = "A" * 28


def _read(biological_clip: str, *, reverse: bool = False, adjacent: int = 20,
          technical_prefix: str = PREFIX) -> pysam.AlignedSegment:
    read = pysam.AlignedSegment()
    read.query_name = "molecule"
    read.flag = 0x10 if reverse else 0
    read.reference_id = 0
    read.reference_start = 100
    read.mapping_quality = 60
    if reverse:
        stored_clip = reverse_complement(biological_clip) + technical_prefix
        read.cigartuples = [(0, adjacent), (4, len(stored_clip))]
        read.query_sequence = "C" * adjacent + stored_clip
    else:
        stored_clip = technical_prefix + biological_clip
        read.cigartuples = [(4, len(stored_clip)), (0, adjacent)]
        read.query_sequence = stored_clip + "C" * adjacent
    return read


def test_full_tso_matches_after_prefix_removal_on_both_orientations() -> None:
    for reverse in (False, True):
        read = _read(TSO, reverse=reverse)
        assert _tso_match_length(
            read, TSO, 8, 2, technical_prefix_bp=28,
            junction_slop_bp=2, min_adjacent_match_bp=6,
        ) == 13


def test_match_crossing_technical_prefix_boundary_is_rejected() -> None:
    # The complete motif exists only across the 28-bp-prefix boundary. The
    # biological portion alone is seven bases and cannot meet min_match=8.
    technical = "C" * 22 + TSO[:6]
    read = _read(TSO[6:], technical_prefix=technical)

    assert not _has_tso(
        read, TSO, 8, 2, technical_prefix_bp=28,
        junction_slop_bp=2, min_adjacent_match_bp=6,
    )


def test_truncated_tso_uses_suffix_and_length_scaled_mismatches() -> None:
    suffix = TSO[-8:]
    one_mismatch = ("A" if suffix[0] != "A" else "C") + suffix[1:]
    two_mismatches = one_mismatch[:-1] + ("A" if one_mismatch[-1] != "A" else "C")

    assert _has_tso(_read(suffix), TSO, 8, 2, technical_prefix_bp=28)
    assert _has_tso(_read(one_mismatch), TSO, 8, 2, technical_prefix_bp=28)
    assert not _has_tso(_read(two_mismatches), TSO, 8, 2, technical_prefix_bp=28)


def test_failed_full_length_clip_cannot_fall_back_to_partial_match() -> None:
    suffix = TSO[-11:]
    mutated_suffix = suffix[:-1] + ("A" if suffix[-1] != "A" else "C")
    # Three total mismatches make the 13-bp observation fail. Discarding the
    # first two bases would leave an 11-bp/one-mismatch suffix, but partial
    # matching is reserved for clips that are genuinely shorter than 13 bp.
    assert not _has_tso(
        _read("AA" + mutated_suffix), TSO, 8, 2,
        technical_prefix_bp=28, junction_slop_bp=2,
    )


def test_tso_must_be_junction_proximal_and_have_an_alignment_anchor() -> None:
    assert not _has_tso(_read(TSO + "AAAAA"), TSO, 8, 2, technical_prefix_bp=28)
    assert not _has_tso(_read(TSO, adjacent=5), TSO, 8, 2, technical_prefix_bp=28)


class _FakeFasta:
    def __init__(self, sequence: str):
        self.sequence = sequence

    def fetch(self, chrom: str, start: int, end: int) -> str:
        return self.sequence[: end - start]


def test_adjacent_genomic_tso_can_be_rejected() -> None:
    read = _read(TSO)

    assert _reference_matches_tso(_FakeFasta(TSO), "chr1", read, TSO, 13, 2)
    assert not _reference_matches_tso(_FakeFasta("C" * 13), "chr1", read, TSO, 13, 2)


def test_adjacent_genomic_check_uses_read_orientation() -> None:
    class OrientationFasta:
        def fetch(self, chrom: str, start: int, end: int) -> str:
            # Reverse reads query downstream of reference_end and normalize the
            # genomic sequence into physical sequencing orientation.
            assert start == 120
            assert end == 133
            return reverse_complement(TSO)

    assert _reference_matches_tso(
        OrientationFasta(), "chr1", _read(TSO, reverse=True), TSO, 13, 2,
    )
