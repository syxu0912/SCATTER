from __future__ import annotations

import sqlite3
from pathlib import Path

from scatter.calling import (
    _adaptive_gap,
    _assign_route_cluster_ids,
    _annotate_reference_tss_support,
    _filter_off_candidates,
    _load_barcode_whitelist,
    _load_utr5,
    _mark_utr5_overlap,
    _merge_mix_candidates,
    _output_paths,
    _route_candidates,
)


def test_output_paths_are_explicit_and_route_aware() -> None:
    mix = _output_paths(Path("/analysis"), "mix")
    assert mix["final_raw_clusters"].name == "mix_route_final_raw_clusters.tsv"
    assert "final_raw_ctss" not in mix
    assert mix["tso_route_raw_clusters"].name == "tso_route_raw_clusters.tsv"
    assert mix["off_route_before_reference_tss_filter_raw_clusters"].name == (
        "off_route_before_reference_tss_filter_raw_clusters.tsv"
    )
    assert mix["off_route_after_reference_tss_filter_raw_clusters"].name == (
        "off_route_after_reference_tss_filter_raw_clusters.tsv"
    )
    assert mix["mix_route_off_after_tso_rescue_raw_clusters"].name == (
        "mix_route_off_after_tso_rescue_raw_clusters.tsv"
    )

    off = _output_paths(Path("/analysis"), "off")
    assert off["final_raw_clusters"].name == "off_route_final_raw_clusters.tsv"
    assert "tso_route_raw_clusters" not in off
    assert "mix_route_off_after_tso_rescue_raw_clusters" not in off


def test_route_cluster_ids_are_deterministic_and_nonempty() -> None:
    rows = [{"gene_id": "G1"}, {"gene_id": "G2"}]

    _assign_route_cluster_ids(rows, "TSO_RAW")

    assert [row["cluster_id"] for row in rows] == [
        "TSO_RAW_00000001", "TSO_RAW_00000002",
    ]


def test_adaptive_gap_is_anchored_at_fifty_reads_and_narrows_with_depth() -> None:
    assert _adaptive_gap(10) == 100
    assert _adaptive_gap(50) == 100
    assert _adaptive_gap(100) == 92
    assert _adaptive_gap(500) == 77
    assert _adaptive_gap(1_000) == 72
    assert _adaptive_gap(10_000) == 59
    assert _adaptive_gap(100_000) == 50
    assert _adaptive_gap(10**30) == 15


def _evidence_db() -> sqlite3.Connection:
    conn = sqlite3.connect(":memory:")
    conn.execute(
        "CREATE TABLE evidence (chrom TEXT, gene TEXT, strand TEXT, pos INTEGER, cb TEXT, umi TEXT, "
        "n_reads INTEGER NOT NULL, tso_reads INTEGER NOT NULL, "
        "PRIMARY KEY(chrom,gene,strand,pos,cb,umi)) WITHOUT ROWID"
    )
    return conn


def test_route_thresholds_are_raw_read_counts() -> None:
    conn = _evidence_db()
    conn.executemany(
        "INSERT INTO evidence VALUES (?,?,?,?,?,?,?,?)",
        [
            ("chr1", "TOO_SMALL", "+", 50, "cellA", "AAAA", 4, 0),
            ("chr1", "PASS", "+", 100, "cellA", "CCCC", 30, 20),
            ("chr1", "PASS", "+", 110, "cellA", "GGGG", 20, 10),
            ("chr1", "PASS", "+", 300, "cellA", "TTTT", 49, 0),
        ],
    )

    rows = _route_candidates(conn, "off", min_gene_reads=5, min_cluster_reads=50)

    assert len(rows) == 1
    assert rows[0]["gene_id"] == "PASS"
    assert rows[0]["start_1"] == 100
    assert rows[0]["end_1"] == 110
    assert rows[0]["qualified_read_count"] == 50
    # The read threshold is intentionally independent of corrected molecules.
    assert rows[0]["corrected_molecule_count"] == 2


def test_tso_route_is_a_subset_of_the_all_read_route() -> None:
    conn = _evidence_db()
    conn.executemany(
        "INSERT INTO evidence VALUES (?,?,?,?,?,?,?,?)",
        [
            ("chr1", "G1", "+", 100, "cellA", "AAAA", 60, 60),
            ("chr1", "G1", "+", 110, "cellB", "CCCC", 70, 0),
        ],
    )

    tso = _route_candidates(conn, "tso", 5, 50)
    off = _route_candidates(conn, "off", 5, 50)

    assert [(row["start_1"], row["qualified_read_count"]) for row in tso] == [(100, 60)]
    assert [(row["start_1"], row["end_1"], row["qualified_read_count"]) for row in off] == [(100, 110, 130)]
    assert tso[0]["tso_read_count"] == 60 and tso[0]["off_read_count"] == 0
    assert off[0]["tso_read_count"] == 60 and off[0]["off_read_count"] == 130


def test_mixed_unions_tso_and_retained_all_read_routes_and_marks_both_sources() -> None:
    conn = _evidence_db()
    conn.executemany(
        "INSERT INTO evidence VALUES (?,?,?,?,?,?,?,?)",
        [
            ("chr1", "G1", "+", 100, "cellA", "AAAA", 60, 60),
            ("chr1", "G1", "+", 110, "cellB", "CCCC", 70, 0),
        ],
    )
    tso = _route_candidates(conn, "tso", 5, 50)
    off = _route_candidates(conn, "off", 5, 50)
    _mark_utr5_overlap(off, {("chr1", "+", "G1"): [(109, 111)]}, 500)
    retained_off = [row for row in off if row["utr5_overlap"]]

    mix = _merge_mix_candidates(tso, retained_off)

    assert len(mix) == 1
    assert mix[0]["start_1"] == 100 and mix[0]["end_1"] == 110
    assert mix[0]["qualified_read_count"] == 130
    assert mix[0]["candidate_source"] == "mix_tso+off"
    assert mix[0]["support_type"] == "both"
    assert mix[0]["tso_read_count"] == 60
    assert mix[0]["off_read_count"] == 130
    assert mix[0]["tso_cluster_count"] == 1
    assert mix[0]["off_cluster_count"] == 1


def test_all_read_cluster_without_utr_overlap_is_not_passed_to_mixed() -> None:
    conn = _evidence_db()
    conn.execute(
        "INSERT INTO evidence VALUES (?,?,?,?,?,?,?,?)",
        ("chr1", "G1", "+", 110, "cellA", "AAAA", 70, 0),
    )
    off = _route_candidates(conn, "off", 5, 50)
    _mark_utr5_overlap(off, {("chr1", "+", "G1"): [(500, 600)]}, 500)

    assert _merge_mix_candidates([], [row for row in off if row["utr5_overlap"]]) == []


def test_utr5_overlap_requires_the_same_gx_gene() -> None:
    rows = [
        {"chrom": "chr1", "strand": "+", "gene_id": "G1", "start_1": 100, "end_1": 110},
        {"chrom": "chr1", "strand": "+", "gene_id": "G2", "start_1": 100, "end_1": 110},
    ]

    _mark_utr5_overlap(rows, {("chr1", "+", "G1"): [(99, 110)]}, 0)

    assert rows[0]["utr5_overlap"] is True
    assert rows[1]["utr5_overlap"] is False


def test_utr5_bed_name_supplies_gene_id_and_coordinates(tmp_path) -> None:
    bed = tmp_path / "utr5.bed"
    bed.write_text("chr1\t100\t120\tENSG000001.7|tx1|5UTR\t0\t+\n")

    assert _load_utr5(bed, 5) == {("chr1", "+", "ENSG000001"): [(95, 125)]}


def test_mixed_keeps_tso_cluster_without_utr5_all_read_support() -> None:
    conn = _evidence_db()
    conn.execute(
        "INSERT INTO evidence VALUES (?,?,?,?,?,?,?,?)",
        ("chr1", "G1", "+", 110, "cellA", "AAAA", 60, 60),
    )
    tso = _route_candidates(conn, "tso", 5, 50)

    mix = _merge_mix_candidates(tso, [])

    assert len(mix) == 1
    assert mix[0]["candidate_source"] == "tso"
    assert mix[0]["support_type"] == "tso"
    assert mix[0]["tso_cluster_count"] == 1
    assert mix[0]["off_cluster_count"] == 0


def test_off_tss_distance_classes_use_summit_and_inclusive_boundaries() -> None:
    rows = [
        {"gene_id": "G1", "chrom": "chr1", "strand": "+", "summit_1": position,
         "start_1": position, "end_1": position}
        for position in (1100, 1500, 2000, 2001)
    ]
    _annotate_reference_tss_support(rows, {("G1", "chr1", "+"): [1000]})

    assert [row["off_tss_class"] for row in rows] == [
        "known_tss", "near_tss_shift", "proximal_atss", "distal_filtered",
    ]
    assert [row["off_tss_distance_pass"] for row in rows] == [True, True, True, False]
    assert [row["nearest_known_tss_distance_bp"] for row in rows] == [100, 500, 1000, 1001]


def test_signed_tss_distance_is_transcription_direction_aware() -> None:
    rows = [
        {"gene_id": "G1", "chrom": "chr1", "strand": "+", "summit_1": 1100,
         "start_1": 1100, "end_1": 1100},
        {"gene_id": "G2", "chrom": "chr1", "strand": "-", "summit_1": 1100,
         "start_1": 1100, "end_1": 1100},
    ]
    reference = {("G1", "chr1", "+"): [1000], ("G2", "chr1", "-"): [1000]}
    _annotate_reference_tss_support(rows, reference)

    assert rows[0]["signed_reference_tss_distance_bp"] == 100
    assert rows[1]["signed_reference_tss_distance_bp"] == -100


def test_distal_off_cluster_is_rescued_only_by_same_locus_tso_overlap() -> None:
    off = [
        {"gene_id": "G1", "chrom": "chr1", "strand": "+", "summit_1": 3000,
         "start_1": 2990, "end_1": 3010, "off_tss_distance_pass": False,
         "off_tss_class": "distal_filtered"},
        {"gene_id": "G1", "chrom": "chr1", "strand": "+", "summit_1": 4000,
         "start_1": 3990, "end_1": 4010, "off_tss_distance_pass": False,
         "off_tss_class": "distal_filtered"},
    ]
    tso = [
        {"gene_id": "G1", "chrom": "chr1", "strand": "+", "summit_1": 3005,
         "start_1": 3005, "end_1": 3005},
        # Coordinate overlap alone is insufficient when the gene differs.
        {"gene_id": "G2", "chrom": "chr1", "strand": "+", "summit_1": 4000,
         "start_1": 3995, "end_1": 4005},
    ]

    retained = _filter_off_candidates(off, tso)

    assert retained == [off[0]]
    assert off[0]["off_route_filter_reason"] == "tso_overlap_rescue"
    assert off[1]["off_route_filter_status"] == "filtered"


def test_all_three_reference_tss_classes_keep_original_threshold() -> None:
    off = [
        {"gene_id": "G1", "chrom": "chr1", "strand": "+", "summit_1": position,
         "start_1": position, "end_1": position}
        for position in (1100, 1500, 2000)
    ]
    _annotate_reference_tss_support(off, {("G1", "chr1", "+"): [1000]})

    retained = _filter_off_candidates(off, [])

    assert retained == off
    assert [row["off_route_filter_reason"] for row in retained] == [
        "known_tss", "near_tss_shift", "proximal_atss",
    ]


def test_standalone_off_stays_tso_independent_while_mixed_can_rescue() -> None:
    row = {
        "gene_id": "G1", "chrom": "chr1", "strand": "+", "summit_1": 3000,
        "start_1": 2990, "end_1": 3010, "off_tss_distance_pass": False,
        "off_tss_class": "distal_filtered",
    }
    tso = [{
        "gene_id": "G1", "chrom": "chr1", "strand": "+", "summit_1": 3005,
        "start_1": 3005, "end_1": 3005,
    }]

    standalone = _filter_off_candidates([dict(row)], [])
    mix = _filter_off_candidates([dict(row)], tso)

    assert standalone == []
    assert len(mix) == 1
    assert mix[0]["off_route_filter_reason"] == "tso_overlap_rescue"


def test_barcode_whitelist_accepts_header_and_first_column(tmp_path) -> None:
    whitelist = tmp_path / "barcodes.tsv"
    whitelist.write_text("cell_id\textra\nAAAC-1\tx\nGGGT-1\ty\n")

    assert _load_barcode_whitelist(whitelist) == {"AAAC-1", "GGGT-1"}
