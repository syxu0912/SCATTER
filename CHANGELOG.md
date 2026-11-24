# Changelog

## 0.5.0

- Added `scatter run` for the complete BAM-to-single-cell TSS/CTSS workflow.
- Packaged and checksum-pinned the hg38 1000 bp sigma=10 Shift model; default cutoff 0.9.
- Replaced per-cluster hard-coded 200 bp scoring with length-aware, batched, deduplicated centered inference.
- Preserved ordered called-cell barcodes and zero-count cells; exported sparse counts, observed usage and 10x-style Matrix Market data.
- Clarified 1-based raw positions versus 0-based half-open AnnData intervals.
- Added paired/R2-only auto profiling, R2-only transferred-TSO checks, gzip text inputs, run provenance and QC report.
- Added executable miniature demonstration and full end-to-end model/count-conservation tests.

Earlier experimental SCATTER 0.4.6 runs remain separate in the research workspace.
