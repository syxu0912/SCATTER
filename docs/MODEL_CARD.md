# Bundled Shift model

| Item | Value |
|---|---|
| ID | `shift_hg38_1000bp_sigma10_n1n2_20261123` |
| Reference genome | human GRCh38/hg38 |
| Input | 1000 bp, ATCG channels |
| Network | optimized sequential CNN + bidirectional LSTM, hidden/center readout, layer normalization |
| Training shift | ten deterministic views per record, Gaussian σ=10 bp; actual retained training views exclude exact heldout matches |
| Negatives | N1 intergenic random; N2 gene-internal hard; N3 excluded |
| Default cutoff | positive-class softmax score ≥0.9 |
| Inference | one centered strand-oriented window, not test-time shift averaging |
| SHA256 | `58410fdd8a913a95f3075dd5c35741b5418d20effab8338d7f329c077c4ce80b` |

For a 1-based representative position `p`, the genomic window is `[p-501,p+499)` in BED coordinates. RNA-minus candidates are reverse complemented. Unsupported characters, missing contigs and windows extending outside the reference are unscorable, not biologically confirmed negatives.

The checkpoint declares its input length, and software dispatches based on that metadata. The default checkpoint is pinned by the packaged JSON manifest. A replacement model is an explicit user choice and must use a compatible checkpoint format and suitable genome/training design.

The cutoff was fixed by the user's analysis policy, not optimized on a new sample's external support labels. Current research comparisons used evidence-supported candidates in K562 and GM12878; these are not independently verified true/false labels. Training/validation sequence, genomic-region and gene exposure affect such comparisons. The demonstration fixture checks software behavior, not predictive performance.

The model scores a reference DNA sequence. It does not infer cell type, transcript abundance, fragment length, or an unsequenced R1 coordinate. UMI count matrices provide the single-cell measurements; model score alone does not establish gene linkage or a precise experimentally verified TSS.
