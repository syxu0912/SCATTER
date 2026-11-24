# Output definitions

## Count matrices

- AnnData: rows are cells, columns are retained TSS clusters or CTSS positions.
- `.X` and `layers['counts']`: sparse corrected-UMI integer counts. UMI directional correction is performed within cell, assigned gene and candidate cluster. Each corrected molecule is assigned to its most-supported CTSS; ties choose the leftmost genomic position.
- TSS and CTSS molecule sums agree. MIX does not add a molecule once for each evidence route.
- Barcode order follows the requested cell list. Cells with zero retained counts remain present, including when no clusters pass.
- Matrix Market: **features × cells**, i.e. transpose of AnnData `.X`. `features.tsv.gz` uses TSS/CTSS IDs in the standard three-column layout; the feature-type token `Gene Expression` is for 10x-reader interoperability and does not mean these features are genes.
- `feature_metadata.tsv.gz` maps feature IDs to genes, positions, model scores and other annotations; `cell_metadata.tsv.gz` supplies cell QC.

## AnnData annotations

| Location | Definition |
|---|---|
| `obs_names` | cell barcode, unchanged |
| `obs.total_counts` | retained corrected UMI count |
| `obs.n_features_by_counts` | retained features with nonzero count |
| `obs.sample_id` | sample label supplied to `scatter run` |
| `var.gene_id` | unique evidence-based GX assignment |
| `var.start_1/end_1/summit_1` | 1-based inclusive cluster bounds / representative position |
| `var.TSS_start/TSS_end`, `Start/End` | 0-based half-open cluster interval |
| `var.Chromosome/Strand` | reference contig and gene/RNA strand |
| `var.deeptssfilter_score` | positive-class model score, not a calibrated biological probability |
| `var.tss_key` | gene/contig/interval/strand key; run-local IDs are not global sample alignments |
| `var.total_counts/n_cells_by_counts` | observed feature UMI total and number of expressing cells |
| `layers.usage` | observed within-cell fractions: by gene for TSS; by cluster for CTSS |
| `uns.scatter` | version, coordinate conventions, model SHA256, input length, threshold and endpoint rule |

`*_two.h5ad` retains all features from genes with at least two features in the corresponding complete matrix; it does not mean exactly two TSSs were chosen per gene. For a method requiring binary TSS pairs, construct explicit pairs before statistical inference.

## Catalogs and provenance

Raw/retained TSVs keep 1-based positions. `sequence_status` distinguishes ambiguous DNA, out-of-bounds windows and missing contigs from valid scored candidates; unscorable candidates do not enter retained matrices. The full raw catalog and metadata remain available.

`run_summary.json` records the concrete model/threshold, evidence layout, TSO source, RNA strand rule, diagnostic limits and result paths. `metadata/input_checksums.json` uses full hashes for the small barcode file and fast size/time/edge fingerprints for large alignment/reference inputs, explicitly named as such. The model itself is always SHA256-verified.

`report.html` and `metadata/matrix_qc.json` summarize cells, TSS features, UMI counts and zero-count cells. `COMPLETE` is written only after successful matrix export and report generation; `metadata/run_status.json` records failed runs with their error.
