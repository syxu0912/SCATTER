#!/usr/bin/env python3
"""Read TSS counts, optionally attach cell metadata and compute a Scanpy embedding."""
import argparse
import anndata as ad
import numpy as np
import pandas as pd

p=argparse.ArgumentParser()
p.add_argument('--input',required=True,help='count/scTSS_count_all.h5ad')
p.add_argument('--output',required=True,help='New analyzed .h5ad path')
p.add_argument('--cell-metadata',help='CSV with cell barcodes as first-column index')
p.add_argument('--embedding',action='store_true',help='Requires scanpy extra; normalize a copy of counts for PCA/UMAP')
p.add_argument('--min-total-umis',type=int,default=1)
a=p.parse_args();x=ad.read_h5ad(a.input)
if a.cell_metadata:
    cells=pd.read_csv(a.cell_metadata,index_col=0)
    if not cells.index.is_unique:raise ValueError('Duplicate metadata cell barcodes')
    if not set(x.obs_names)&set(cells.index):raise ValueError('No matching cell barcodes in metadata')
    for c in cells.columns:
        if c in x.obs:raise ValueError(f'Column already exists: {c}')
        x.obs[c]=cells[c].reindex(x.obs_names)
if a.embedding:
    import scanpy as sc
    sc.pp.filter_cells(x,min_counts=a.min_total_umis)
    sc.pp.filter_genes(x,min_cells=1)
    if min(x.shape)<3:raise ValueError('Too few nonzero cells/features for an embedding; read counts without --embedding')
    sc.pp.normalize_total(x,target_sum=1e4);sc.pp.log1p(x)
    sc.pp.highly_variable_genes(x,n_top_genes=min(2000,x.n_vars))
    features=int(x.var.highly_variable.sum());components=min(30,x.n_obs-1,features-1)
    if components<1:raise ValueError('Too few variable features')
    sc.tl.pca(x,n_comps=components,svd_solver='arpack',use_highly_variable=True)
    sc.pp.neighbors(x,n_neighbors=min(15,x.n_obs-1),n_pcs=components)
    sc.tl.umap(x)
x.write_h5ad(a.output,compression='gzip')
print(x)
