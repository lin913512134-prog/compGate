# Input and output schema

## Expression matrix

A genes × samples table. In Python, use a `pandas.DataFrame` with gene symbols as index and sample IDs as columns. CSV readers should use the first column as the gene index. Values must be finite and on a consistent expression scale.

## Group vector

A binary numeric array aligned to expression columns: case = 1 and control = 0. The original case/control order must be preserved; do not silently reorder samples.

## Fraction matrix

A samples × cell-types `DataFrame` aligned to expression columns. Fractions must be finite and sum to one within tolerance. The joint model removes zero-variance fractions, drops the final active fraction as reference and reports additional collinearity handling, fitted rank and residual degrees of freedom.

## Signature registry

Required fields: `signature`, `genes`, `list_status` and `counts_toward_30`. Multiple genes are pipe-delimited. At least two measured genes are required for the paired workflows. Source citations are required for the public 30-signature registry.

## Main quickstart outputs

- `quickstart.json`: package version, synthetic scope, attribution record and joint-model state.
- `local_verification.json`: test counts, quickstart failures, Figure 4 anchors and environment.

## Full frozen-data entry points

- `R18_fraction_error_per_run.csv` → benchmark scenario×method medians/means/minima/maxima and 36-row factorial coverage.
- `R18_reference_convention_long.csv` → paired per-reference cells, per-signature stalls, reference summary and 50 common six-arm cells.
- `examples/figure4_source/` → Figure 4 PNG/PDF and `figure4_anchors.json`.

## Public roles

The source repository contains executable code and minimal figure inputs. The Zenodo candidate contains the manuscript, five figures, source data, two main tables, five supplementary entry points, verification logs and reviewer manifest.
