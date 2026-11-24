from __future__ import annotations

import csv
from pathlib import Path
import sqlite3

import anndata as ad

from scatter.calling import _directional_umi_correction
from scatter.matrix import export_filtered_h5ad


CLUSTER_FIELDS = [
    "cluster_id", "gene_id", "chrom", "start_1", "end_1", "summit_1", "strand",
    "corrected_molecule_count", "deeptssfilter_score", "candidate_source", "filter_status",
]


def _write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_export_filtered_h5ad(tmp_path: Path) -> None:
    database = tmp_path / "evidence.sqlite"
    conn = sqlite3.connect(database)
    conn.execute(
        "CREATE TABLE evidence (chrom TEXT, gene TEXT, strand TEXT, pos INTEGER, cb TEXT, umi TEXT, "
        "n_reads INTEGER NOT NULL, tso_reads INTEGER NOT NULL, "
        "PRIMARY KEY(chrom,gene,strand,pos,cb,umi)) WITHOUT ROWID"
    )
    conn.executemany(
        "INSERT INTO evidence VALUES (?,?,?,?,?,?,?,?)",
        [
            ("chr1", "G1", "+", 100, "cellA", "AAAA", 4, 0),
            ("chr1", "G1", "+", 101, "cellA", "AAAT", 1, 0),
            ("chr1", "G1", "+", 101, "cellA", "CCCC", 1, 0),
            ("chr1", "G1", "+", 100, "cellB", "GGGG", 1, 0),
            ("chr1", "G1", "+", 200, "cellA", "TTTT", 1, 0),
            ("chr1", "G2", "+", 300, "cellC", "ACGT", 5, 1),
        ],
    )
    conn.commit()
    conn.close()

    clusters = tmp_path / "high_confidence_clusters.tsv"
    _write_tsv(
        clusters,
        CLUSTER_FIELDS,
        [
            {"cluster_id": "C1", "gene_id": "G1", "chrom": "chr1", "start_1": 100, "end_1": 101,
             "summit_1": 100, "strand": "+", "corrected_molecule_count": 3,
             "deeptssfilter_score": 0.99, "candidate_source": "off", "filter_status": "pass"},
            {"cluster_id": "C2", "gene_id": "G1", "chrom": "chr1", "start_1": 200, "end_1": 200,
             "summit_1": 200, "strand": "+", "corrected_molecule_count": 1,
             "deeptssfilter_score": 0.98, "candidate_source": "off", "filter_status": "pass"},
            {"cluster_id": "C3", "gene_id": "G2", "chrom": "chr1", "start_1": 300, "end_1": 300,
             "summit_1": 300, "strand": "+", "corrected_molecule_count": 1,
             "deeptssfilter_score": 0.97, "candidate_source": "tso", "filter_status": "pass"},
        ],
    )
    genes = tmp_path / "genes.tsv"
    _write_tsv(
        genes,
        ["gene_id", "gene_name", "gene_type", "chrom", "start_1", "end_1", "strand"],
        [
            {"gene_id": "G1", "gene_name": "Gene1", "gene_type": "protein_coding", "chrom": "chr1",
             "start_1": 1, "end_1": 250, "strand": "+"},
            {"gene_id": "G2", "gene_name": "Gene2", "gene_type": "lncRNA", "chrom": "chr1",
             "start_1": 251, "end_1": 400, "strand": "+"},
        ],
    )

    result = export_filtered_h5ad(
        evidence_db=database,
        clusters_tsv=clusters,
        genes_tsv=genes,
        outdir=tmp_path / "out",
        umi_corrector=_directional_umi_correction,
    )

    assert result["tss_cells"] == 3
    assert result["tss_features"] == 3
    assert result["ctss_features"] == 4
    tss = ad.read_h5ad(tmp_path / "out" / "count" / "scTSS_count_all.h5ad")
    assert tss.X.toarray().sum() == 5
    assert tss.var.loc["C1", "gene_name"] == "Gene1"
    assert tss[tss.obs_names == "cellA", tss.var_names == "C1"].X[0, 0] == 2
    tss_two = ad.read_h5ad(tmp_path / "out" / "count" / "scTSS_count_two.h5ad")
    assert tss_two.var_names.tolist() == ["C1", "C2"]
    ctss_two = ad.read_h5ad(tmp_path / "out" / "ctss" / "all_ctss_two.h5ad")
    assert ctss_two.n_vars == 3


def test_export_empty_filtered_set(tmp_path: Path) -> None:
    database = tmp_path / "evidence.sqlite"
    conn = sqlite3.connect(database)
    conn.execute(
        "CREATE TABLE evidence (chrom TEXT, gene TEXT, strand TEXT, pos INTEGER, cb TEXT, umi TEXT, "
        "n_reads INTEGER NOT NULL, tso_reads INTEGER NOT NULL, "
        "PRIMARY KEY(chrom,gene,strand,pos,cb,umi)) WITHOUT ROWID"
    )
    conn.close()
    clusters = tmp_path / "clusters.tsv"
    genes = tmp_path / "genes.tsv"
    _write_tsv(clusters, CLUSTER_FIELDS, [])
    _write_tsv(genes, ["gene_id", "gene_name", "gene_type", "chrom", "start_1", "end_1", "strand"], [])

    export_filtered_h5ad(
        evidence_db=database,
        clusters_tsv=clusters,
        genes_tsv=genes,
        outdir=tmp_path / "out",
        umi_corrector=_directional_umi_correction,
    )
    empty = ad.read_h5ad(tmp_path / "out" / "count" / "scTSS_count_all.h5ad")
    assert empty.shape == (0, 0)
    assert "gene_id" in empty.var.columns
