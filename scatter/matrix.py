"""Export filtered SCATTER calls as CamoTSS-compatible AnnData matrices."""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
import gzip
from pathlib import Path
import sqlite3
from typing import Callable

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

from . import __version__
from .io import read_barcodes


UmiCorrector = Callable[[dict[str, int]], dict[str, str]]

_INTEGER_CLUSTER_COLUMNS = {
    "start_1", "end_1", "summit_1", "adaptive_gap_bp", "qualified_read_count", "raw_site_umi_count",
    "corrected_molecule_count", "ctss_count", "nearest_known_tss_distance_bp",
    "tso_matched_alignment_count", "tso_total_alignment_count", "tso_read_count", "off_read_count",
    "tso_cluster_count", "off_cluster_count", "source_cluster_count", "utr5_flank_bp",
    "nearest_reference_tss_1", "signed_reference_tss_distance_bp",
    "parent_start_1", "parent_end_1", "parent_qualified_read_count", "cluster_max_width_bp",
    "model_sequence_length_bp",
}
_FLOAT_CLUSTER_COLUMNS = {"deeptssfilter_score", "tso_fraction_of_alignments"}


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _load_gene_metadata(path: Path) -> dict[str, dict[str, str]]:
    return {row["gene_id"].split(".")[0]: row for row in _read_tsv(path)}


def _selected_read_count(row: sqlite3.Row, candidate_source: str) -> int:
    # A mix_tso+off candidate is represented by its all-read (off-route)
    # member in calling._merge_mix_candidates. TSO-only candidates retain
    # candidate_source=tso.
    return int(row["tso_reads"] if candidate_source == "tso" else row["n_reads"])


def _quantify_cluster(
    conn: sqlite3.Connection,
    cluster: dict[str, str],
    umi_corrector: UmiCorrector,
) -> tuple[Counter[str], dict[int, Counter[str]]]:
    rows = conn.execute(
        "SELECT cb,umi,pos,n_reads,tso_reads FROM evidence "
        "WHERE chrom=? AND gene=? AND strand=? AND pos BETWEEN ? AND ?",
        (
            cluster["chrom"], cluster["gene_id"], cluster["strand"],
            int(cluster["start_1"]), int(cluster["end_1"]),
        ),
    ).fetchall()
    per_cell: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    per_position: dict[tuple[str, str], Counter[int]] = defaultdict(Counter)
    source = cluster.get("candidate_source", "off")
    for row in rows:
        selected = _selected_read_count(row, source)
        if not selected:
            continue
        cb, umi, position = str(row["cb"]), str(row["umi"]), int(row["pos"])
        per_cell[cb][umi] += selected
        per_position[(cb, umi)][position] += selected

    cluster_counts: Counter[str] = Counter()
    ctss_counts: dict[int, Counter[str]] = defaultdict(Counter)
    for cb, abundances in per_cell.items():
        parents = umi_corrector(abundances)
        corrected_positions: dict[str, Counter[int]] = defaultdict(Counter)
        for umi, parent in parents.items():
            corrected_positions[parent].update(per_position[(cb, umi)])
        for counts in corrected_positions.values():
            # Match the CTSS assignment used by the caller: most supporting
            # alignments first, then the left-most genomic position on ties.
            position = min(counts, key=lambda value: (-counts[value], value))
            cluster_counts[cb] += 1
            ctss_counts[position][cb] += 1
    return cluster_counts, dict(ctss_counts)


def _coerce_cluster_value(column: str, value: str | None) -> object:
    if value in (None, ""):
        return np.nan if column in _INTEGER_CLUSTER_COLUMNS | _FLOAT_CLUSTER_COLUMNS else ""
    if column in _INTEGER_CLUSTER_COLUMNS:
        return int(value)
    if column in _FLOAT_CLUSTER_COLUMNS:
        return float(value)
    return value


def _cluster_var(
    clusters: list[dict[str, str]],
    gene_metadata: dict[str, dict[str, str]],
) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for cluster in clusters:
        record = {key: _coerce_cluster_value(key, value) for key, value in cluster.items() if key != "cluster_id"}
        gene = gene_metadata.get(cluster["gene_id"], {})
        # These aliases mirror the fields most commonly consumed from a
        # CamoTSS scTSS_count_all.h5ad file.
        record.update(
            {
                "TSS_start": int(cluster["start_1"]) - 1,
                "TSS_end": int(cluster["end_1"]),
                "Chromosome": cluster["chrom"],
                "Strand": cluster["strand"],
                "gene_name": gene.get("gene_name", ""),
                "gene_type": gene.get("gene_type", ""),
                "Feature": "TSS",
                "Start": int(cluster["start_1"]) - 1,
                "End": int(cluster["end_1"]),
                "len": int(cluster["end_1"]) - int(cluster["start_1"]) + 1,
                "tss_key": f"{cluster['gene_id']}|{cluster['chrom']}:{cluster['start_1']}-{cluster['end_1']}:{cluster['strand']}",
                "gene_start_1": int(gene["start_1"]) if gene.get("start_1") else np.nan,
                "gene_end_1": int(gene["end_1"]) if gene.get("end_1") else np.nan,
            }
        )
        records.append(record)
    index = pd.Index([row["cluster_id"] for row in clusters], name="cluster_id")
    if not records:
        return pd.DataFrame(
            index=index,
            columns=[
                "gene_id", "chrom", "start_1", "end_1", "summit_1", "strand",
                "deeptssfilter_score", "filter_status", "TSS_start", "TSS_end",
                "Chromosome", "Strand", "gene_name", "gene_type", "gene_start_1", "gene_end_1",
            ],
        )
    return pd.DataFrame(records, index=index)


def _ctss_var(
    features: list[tuple[dict[str, str], int, int]],
    gene_metadata: dict[str, dict[str, str]],
) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    names: list[str] = []
    for cluster, position, molecule_count in features:
        gene = gene_metadata.get(cluster["gene_id"], {})
        names.append(f"{cluster['cluster_id']}|{cluster['chrom']}:{position}:{cluster['strand']}")
        records.append(
            {
                "cluster_id": cluster["cluster_id"],
                "gene_id": cluster["gene_id"],
                "Chromosome": cluster["chrom"],
                "CTSS": position,
                "Strand": cluster["strand"],
                "corrected_molecule_count": molecule_count,
                "gene_name": gene.get("gene_name", ""),
                "gene_type": gene.get("gene_type", ""),
            }
        )
    index = pd.Index(names, name="ctss_id")
    if not records:
        return pd.DataFrame(
            index=index,
            columns=[
                "cluster_id", "gene_id", "Chromosome", "CTSS", "Strand",
                "corrected_molecule_count", "gene_name", "gene_type",
            ],
        )
    return pd.DataFrame(records, index=index)


def _make_adata(
    feature_names: list[str],
    counts: dict[tuple[str, str], int],
    var: pd.DataFrame,
    matrix_type: str,
    cells: list[str] | None = None,
) -> ad.AnnData:
    cells = cells if cells is not None else sorted({cell for cell, _ in counts})
    if len(cells) != len(set(cells)):
        raise ValueError('Duplicate matrix cell barcodes')
    if not {cell for cell, _ in counts}.issubset(cells):
        raise ValueError('Evidence contains cells outside the requested barcode list')
    cell_index = {cell: index for index, cell in enumerate(cells)}
    feature_index = {feature: index for index, feature in enumerate(feature_names)}
    coordinates = [
        (cell_index[cell], feature_index[feature], count)
        for (cell, feature), count in counts.items()
    ]
    if coordinates:
        row_indices, column_indices, values = zip(*coordinates)
        matrix = sparse.csr_matrix(
            (np.asarray(values, dtype=np.int32), (row_indices, column_indices)),
            shape=(len(cells), len(feature_names)),
            dtype=np.int32,
        )
    else:
        matrix = sparse.csr_matrix((len(cells), len(feature_names)), dtype=np.int32)
    result = ad.AnnData(
        X=matrix,
        obs=pd.DataFrame(index=pd.Index(cells, name="cell_barcode")),
        var=var,
    )
    result.uns["scatter"] = {
        "version": __version__,
        "matrix_type": matrix_type,
        "scope": "DeepTSSFilter-passing, uniquely assigned clusters",
        "count_unit": "directionally corrected UMI molecules",
        "coordinates": "start_1/end_1/summit_1/CTSS are 1-based; TSS_start/TSS_end and Start/End are 0-based half-open",
    }
    result.layers['counts'] = matrix.copy()
    result.obs['total_counts'] = np.asarray(matrix.sum(axis=1)).ravel()
    result.obs['n_features_by_counts'] = matrix.getnnz(axis=1)
    result.var['total_counts'] = np.asarray(matrix.sum(axis=0)).ravel()
    result.var['n_cells_by_counts'] = matrix.getnnz(axis=0)
    grouping = 'gene_id' if matrix_type == 'TSS cluster' else 'cluster_id'
    result.layers['usage'] = within_group_usage(matrix, var[grouping].tolist())
    result.uns['scatter']['usage_group'] = grouping
    result.uns['scatter']['usage_definition'] = 'Within-cell fractions among retained features; not inferred differential usage or missing-data imputation'
    return result


def within_group_usage(matrix, groups):
    """Sparse within-cell feature usage, without densifying cells x features."""
    matrix = matrix.tocsr()
    if matrix.shape[1] == 0 or matrix.nnz == 0:
        return sparse.csr_matrix(matrix.shape, dtype=np.float32)
    codes, unique = pd.factorize(groups, sort=False)
    mapping = sparse.csr_matrix((np.ones(len(codes), dtype=np.int32), (np.arange(len(codes)), codes)),
                                 shape=(len(codes), len(unique)))
    totals = (matrix @ mapping).tocsr()
    coo = matrix.tocoo()
    denominator = np.asarray(totals[coo.row, codes[coo.col]]).ravel()
    if (denominator <= 0).any():
        raise ValueError('Invalid within-group count denominator')
    return sparse.csr_matrix((coo.data.astype(np.float32)/denominator, (coo.row, coo.col)), shape=matrix.shape)


def write_mtx(adata, directory):
    """10x/Seurat orientation: features x cells; AnnData remains cells x features."""
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    with gzip.open(directory/'matrix.mtx.gz', 'wb') as handle:
        mmwrite(handle, adata.X.T.tocoo(), field='integer', symmetry='general')
    with gzip.open(directory/'barcodes.tsv.gz', 'wt') as handle:
        handle.write(''.join(str(x)+'\n' for x in adata.obs_names))
    with gzip.open(directory/'features.tsv.gz', 'wt') as handle:
        for identifier in adata.var_names:
            # Gene Expression is the interoperable 10x feature-type label.
            handle.write(f'{identifier}\t{identifier}\tGene Expression\n')
    adata.var.to_csv(directory/'feature_metadata.tsv.gz', sep='\t', compression='gzip')
    adata.obs.to_csv(directory/'cell_metadata.tsv.gz', sep='\t', compression='gzip')


def _write_pair(adata: ad.AnnData, all_path: Path, two_path: Path) -> None:
    all_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(all_path, compression="gzip")
    gene_counts = adata.var["gene_id"].value_counts()
    genes_with_multiple_features = set(gene_counts[gene_counts >= 2].index)
    subset = adata[:, adata.var["gene_id"].isin(genes_with_multiple_features)].copy()
    subset.write_h5ad(two_path, compression="gzip")


def export_filtered_h5ad(
    *,
    evidence_db: str | Path,
    clusters_tsv: str | Path,
    genes_tsv: str | Path,
    outdir: str | Path,
    umi_corrector: UmiCorrector,
    barcode_whitelist: str | Path | None = None,
    sample_id: str | None = None,
    run_metadata: dict | None = None,
) -> dict[str, object]:
    """Write filtered TSS and CTSS matrices from an existing SCATTER run."""
    evidence_path = Path(evidence_db).expanduser().resolve()
    cluster_path = Path(clusters_tsv).expanduser().resolve()
    gene_path = Path(genes_tsv).expanduser().resolve()
    root = Path(outdir).expanduser().resolve()
    for required in (evidence_path, cluster_path, gene_path):
        if not required.is_file():
            raise FileNotFoundError(f"Required matrix-export input does not exist: {required}")

    clusters = [row for row in _read_tsv(cluster_path) if row.get("filter_status", "pass") == "pass"]
    gene_metadata = _load_gene_metadata(gene_path)
    tss_counts: dict[tuple[str, str], int] = {}
    ctss_counts: dict[tuple[str, str], int] = {}
    ctss_features: list[tuple[dict[str, str], int, int]] = []

    conn = sqlite3.connect(f"file:{evidence_path}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    try:
        cells = (read_barcodes(barcode_whitelist) if barcode_whitelist is not None else
                 [row[0] for row in conn.execute('SELECT DISTINCT cb FROM evidence ORDER BY cb')])
        for cluster in clusters:
            per_cluster, per_ctss = _quantify_cluster(conn, cluster, umi_corrector)
            observed_total = sum(per_cluster.values())
            expected_total = int(cluster["corrected_molecule_count"])
            if observed_total != expected_total:
                raise ValueError(
                    f"Matrix count mismatch for {cluster['cluster_id']}: "
                    f"reconstructed {observed_total}, cluster table reports {expected_total}"
                )
            for cell, count in per_cluster.items():
                tss_counts[(cell, cluster["cluster_id"])] = count
            for position in sorted(per_ctss):
                feature_name = f"{cluster['cluster_id']}|{cluster['chrom']}:{position}:{cluster['strand']}"
                molecule_count = sum(per_ctss[position].values())
                ctss_features.append((cluster, position, molecule_count))
                for cell, count in per_ctss[position].items():
                    ctss_counts[(cell, feature_name)] = count
    finally:
        conn.close()

    tss_feature_names = [cluster["cluster_id"] for cluster in clusters]
    ctss_feature_names = [
        f"{cluster['cluster_id']}|{cluster['chrom']}:{position}:{cluster['strand']}"
        for cluster, position, _ in ctss_features
    ]
    tss = _make_adata(tss_feature_names, tss_counts, _cluster_var(clusters, gene_metadata), "TSS cluster", cells)
    ctss = _make_adata(ctss_feature_names, ctss_counts, _ctss_var(ctss_features, gene_metadata), "CTSS", cells)
    for matrix in (tss,ctss):
        if sample_id is not None: matrix.obs['sample_id']=sample_id
        matrix.uns['scatter'].update(run_metadata or {})

    tss_all = root / "count" / "scTSS_count_all.h5ad"
    tss_two = root / "count" / "scTSS_count_two.h5ad"
    ctss_all = root / "ctss" / "all_ctss.h5ad"
    ctss_two = root / "ctss" / "all_ctss_two.h5ad"
    _write_pair(tss, tss_all, tss_two)
    _write_pair(ctss, ctss_all, ctss_two)
    write_mtx(tss, root/'count'/'matrix')
    write_mtx(ctss, root/'ctss'/'matrix')
    return {
        "tss_cells": tss.n_obs,
        "tss_features": tss.n_vars,
        "ctss_cells": ctss.n_obs,
        "ctss_features": ctss.n_vars,
        "files": {
            "scTSS_count_all": str(tss_all),
            "scTSS_count_two": str(tss_two),
            "all_ctss": str(ctss_all),
            "all_ctss_two": str(ctss_two),
        },
    }
