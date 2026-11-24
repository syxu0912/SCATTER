# SCATTER

SCATTER discovers transcription-start-site (TSS) clusters from tagged 5′ single-cell BAM files, filters candidates with a bundled DNA-sequence model, and counts UMI molecules for each cell and TSS/CTSS. It supports the paired and R2-only BAM conventions used in this project.

**Default filter:** hg38/GRCh38, 1000 bp Shift model, training σ=10 bp, N1/N2 negatives, score ≥0.9. Inference uses one centered window. The weight file is included in the Python distribution and verified by SHA256.

## Install

Use Python 3.10 or newer in an isolated environment. For CPU inference, install a suitable [PyTorch CPU distribution](https://pytorch.org/get-started/locally/) first if needed, then:

```bash
python -m pip install .
scatter doctor
```

The repository is private; authenticate with GitHub before cloning:

```bash
git clone https://github.com/syxu0912/SCATTER.git
cd SCATTER
python -m pip install .
```

No research genome, BAM, credentials, or experiment output is included in the repository; only a miniature demonstration dataset is bundled.

## Try the full workflow

```bash
scatter demo --outdir demo-paired --layout paired --mode mix
scatter demo --outdir demo-r2 --layout r2only --mode off
```

Both commands run the real bundled model and produce usable cell-by-TSS and cell-by-CTSS matrices. The demo uses artificial cells/read alignments and a small chromosome containing GRCh38-derived positive DNA windows. It is a software/count-conservation check, not an accuracy benchmark.

## Run on your BAM

Use the FASTA and GTF matching the BAM alignment reference. Supply the **called-cell** barcode list, such as `filtered_feature_bc_matrix/barcodes.tsv.gz`.

```bash
scatter run \
  --bam sample.bam \
  --fasta genome.fa \
  --gtf genes.gtf.gz \
  --barcodes barcodes.tsv.gz \
  --layout paired \
  --mode mix \
  --outdir results/sample \
  --threads 4 --device auto
```

For the corrected R2-only endpoint convention:

```bash
scatter run --bam r2only.bam --fasta genome.fa --gtf genes.gtf \
  --barcodes barcodes.tsv.gz --layout r2only --mode off --outdir results/r2
```

`--layout auto` infers the mate layout from alignment flags; explicit `paired` or `r2only` is preferable for reproducible runs. `--all-barcodes` is an explicit alternative to a called-cell list; it does **not** perform cell calling.

### Three evidence modes

| Mode | Candidate evidence and selection |
|---|---|
| `tso` | TSO-positive evidence; no reference-TSS distance restriction |
| `off` | All qualifying reads, without TSO filtering; retain candidates within the same gene's reference TSS ±1000 bp |
| `mix` | Independently call TSO and all-read routes, allow overlapping TSO evidence to rescue off-route candidates, and merge overlapping candidates without counting the same evidence twice |

The modes select read evidence; **all modes output both TSS-cluster and CTSS count matrices**. Defaults require 50 qualifying reads per route/gene and 50 per raw cluster, followed by Shift filtering. These are read thresholds, not UMI thresholds; the reported matrices contain directionally corrected UMI molecule counts. Cluster width is capped at 200 bp using the existing valley-aware splitting rule.

### BAM contract

- Coordinate-sorted BAM; mapped primary records, MAPQ ≥30 by default.
- Corrected `CB` and `UB` tags, and one unambiguous `GX` gene assignment matching the GTF. Raw `CR`/`UR` tags can be selected explicitly with a whitelist; SCATTER does not correct cell barcodes itself.
- Paired mode uses the physical R1 5′ endpoint. The existing long-R1 protocol check requires the dominant R1 length to be at least 100 bp.
- R2-only mode stores the **gene/RNA strand**, using the alignment's left boundary for RNA `+`, right boundary for RNA `−`. It does not infer an absent R1 start from fragment length. A paired-flag BAM can also select R2 explicitly.
- For `tso`/`mix`, `--tso-source auto` prefers a complete boolean `ZT` tag; otherwise paired R1 uses the physical 5′ soft clip. If soft clips retain a technical prefix, set `--technical-prefix-bp` accordingly.
- **R2-only `tso`/`mix` requires transferred TSO evidence**, normally `ZT`. R2 alone cannot recover a TSO removed with R1. Missing transferred evidence is rejected rather than silently running another mode.

BAM/FASTA indices are reused when present; newly needed indices are written inside the output directory. Source BAM/FASTA/GTF files are not changed. A nonempty output directory is not overwritten.

## Outputs

```text
results/sample/
  count/scTSS_count_all.h5ad       cells × retained TSS clusters
  count/scTSS_count_two.h5ad       genes with at least two retained TSS clusters
  count/matrix/                   feature × cell Matrix Market + barcodes/features
  ctss/all_ctss.h5ad              cells × CTSS features
  ctss/all_ctss_two.h5ad           genes with at least two retained CTSS features
  ctss/matrix/                    CTSS Matrix Market files
  clusters/                      raw, route-specific, and model-passing catalogs
  metadata/                      configuration, input fingerprints, QC, model/endpoint provenance
  report.html                    human-readable run summary
  run_summary.json               machine-readable summary
  work/ctss_evidence.sqlite       per-cell UMI evidence for reproducible recounting
  COMPLETE                       written only after successful completion
```

AnnData `.X` and `.layers['counts']` are sparse **integer corrected-UMI counts**. `.layers['usage']` contains observed within-cell fractions: within gene for TSS clusters, within cluster for CTSS. This layer is not imputed or a differential-usage test. `.obs` has per-cell count QC; `.var` contains gene, coordinates, route and model score. Requested barcode order and zero-count cells are retained in every matrix.

Raw table `start_1/end_1/summit_1` and CTSS coordinates are 1-based. AnnData `TSS_start/TSS_end`, `Start/End`, and BED-style intervals are **0-based half-open**. See [output definitions](docs/OUTPUTS.md).

```python
import anndata as ad
a = ad.read_h5ad('results/sample/count/scTSS_count_all.h5ad')
print(a.shape, a.obs.head(), a.var[['gene_id', 'summit_1', 'deeptssfilter_score']].head())
```

See [中文使用说明](docs/README_CN.md), [downstream examples](examples/analyze_tss.py), [model card](docs/MODEL_CARD.md), and `scatter run --help`.

## Interpretation and scope

The model scores reference DNA sequence, not cell-specific expression. Identical chromosome/position/strand loci therefore have identical model scores across cells or samples. A model-passing cluster is a candidate supported by the trained sequence model, not an experimentally confirmed TSS.

The default model was trained on human hg38; other species/assemblies require an appropriately trained replacement checkpoint. TSO/OFF/MIX, UMI correction and model scoring are implemented by SCATTER; CamoTSS is a design/reference point for the single-cell output workflow, not a required executable dependency. CAGE/GRO-cap/RAMPAGE/PLS files used in the research evaluation are not required to run the software and are not used to gate the default model filter.

Independent per-sample cluster IDs are local to each run. Do not concatenate separate matrices by `SCATTER_...` IDs as if they were a shared TSS catalog. Use genomic/gene annotations and a common re-quantified catalog for cross-sample differential TSS analysis.

## Development

```bash
python -m pip install -e '.[test]'
python -m pytest -q
python -m build
```

Tests cover read selection, endpoints/strands, TSO matching, route merging, cluster splitting, UMI conservation, sparse usage, coordinates, zero-count cells and all six paired/R2-only × mode combinations with the actual packaged model.

## Reference and license

The output organization is inspired by [CamoTSS](https://github.com/StatBiomed/CamoTSS) and its [documentation](https://camotss.readthedocs.io/en/latest/run_CamoTSS.html). Cite the relevant tools and methods used in your analysis. SCATTER is distributed under Apache-2.0; see `LICENSE` and `NOTICE`.
