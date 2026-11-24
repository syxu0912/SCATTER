# Release 0.5.0 validation

Local validation on Python 3.12 passed **71 tests**. These include six complete runs through the actual bundled Shift model: paired/R2-only BAM × `tso`/`off`/`mix`.

The miniature example has three requested cells, including one with zero retained counts. Expected conserved corrected-UMI counts are:

| Mode | Retained TSS features | Total TSS UMIs | Total CTSS UMIs |
|---|---:|---:|---:|
| `off` | 3 | 7 | 7 |
| `mix` | 3 | 7 | 7 |
| `tso` | 2 | 5 | 5 |

The same results hold for both input layouts. Tests verify the whitelist cell order, zero-count cell, directional UMI correction, antisense gene assignment, source BAM/FASTA immutability, RNA-facing R2-only endpoints, sparse usage, integer matrix orientation and model checksum. gzip GTF/barcode input, invalid sequence windows and missing transferred TSO tags are also checked.

The built wheel was installed into an isolated environment outside the source tree. `scatter doctor`, paired/MIX and R2-only/OFF demonstrations passed from that installed package. This checks that the model and example assets are actually present in the distribution.

Bounded integrations using real paired and R2-only BAM subsets were also checked locally, including TSS/CTSS UMI conservation and identical cell indexing. Research BAMs, barcode lists and resulting matrices are not included in this repository.

These tests establish software/count consistency, not predictive accuracy or an independent biological benchmark. The default model's research training/validation exposure is relevant when evaluating real-sample predictions.

GitHub Actions runs the test suite, builds the distribution and runs an installed demonstration on Python 3.10 and 3.12. See the repository's Actions tab for current remote results; the local report above does not substitute for a pending remote job.
