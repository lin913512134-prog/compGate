# Changelog

## 0.4.0 — 2026-11-02

- Added `compgate.joint_model` using the exact pure-function implementation used for the revised CMPB analyses.
- Added portable joint-model tests, deterministic quickstart and Figure 4 reconstruction.
- Added frozen-data benchmark and paired-summary entry points.
- Added locked reproduction dependencies, input/output schemas and independent-reproduction record form.
- Preserved the compGate 0.3.6 scoring and attribution-report implementation.
- Public source-repository release; PyPI, a new Zenodo DOI, Linux execution and independent-team validation are not claimed by this entry.

## 0.4.0 — 2026-11-06 (asset-only correction)

- Figures 1, 2, 3 and 5 were rebuilt for the editorial revision; Figure 4 is unchanged.
- `scripts/rebuild_all_figures.py` is now the single self-contained source that reproduces all five figures from the frozen `source_data/` CSVs. It writes no Source Data and has no absolute-path dependency.
- Figure assets are distributed in two formats only: **vector master = PDF; preview = PNG**. EPS and TIFF were retired as redundant.
- `reviewer_manifest.csv` and `SHA256SUMS.txt` were recomputed; the `PENDING_NEW_ZENODO_VERSION` placeholders were eliminated (18 resolved to the DOI, 7 rows retired with EPS/TIFF) `10.5281/zenodo.23889776`.
- Python package code is unchanged; the PyPI 0.4.0 wheel and sdist are unaffected and no version bump was made.
- The frozen `verification/quickstart.json` and `verification/local_verification.json` are retained verbatim as the 2026-11-03 Python 3.12.15 re-run record. Their embedded figure PNG hashes describe the pre-revision figures; the current figure hashes are recorded in `verification/all_figures/figure_build_ledger.json`. See the PD5 section of the Provenance supplement for the two quickstart cost measurements.


## 0.3.6

- Deposited scoring and `attribution_report()` implementation; not published to PyPI as 0.3.6.

## 0.3.5

- Last public PyPI version verified on 2026-11-02.
