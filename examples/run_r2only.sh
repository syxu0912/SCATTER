#!/usr/bin/env bash
set -euo pipefail
: "${BAM:?Set BAM}" "${FASTA:?Set FASTA}" "${GTF:?Set GTF}" "${BARCODES:?Set BARCODES}" "${OUTDIR:?Set OUTDIR}"
scatter run --bam "$BAM" --fasta "$FASTA" --gtf "$GTF" --barcodes "$BARCODES" \
  --layout r2only --mode off --outdir "$OUTDIR" --device auto
