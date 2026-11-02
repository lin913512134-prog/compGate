# compGate 0.4.0

Portable software for **A joint-model framework and benchmark for composition attribution in blood transcriptomics**.
It combines the existing compGate scoring/reporting API with `compgate.joint_model`, the joint-model implementation used for the revised CMPB analyses.

**Release status:** version 0.4.0 was prepared for the CMPB source release on 2026-11-02. This repository provides the public source release. PyPI and Zenodo publication, a version-specific DOI, Linux execution and independent-team validation are separate checks and are not claimed by this source release.

## Included functionality

- `compgate.recap_gate`, `gate_report`, `attribution_report`, `CAR` and `mean_z`.
- `compgate.joint_model`: common-scale joint diagnostics, FWL checks, rank/degrees-of-freedom states, dual scoring, training-fixed transformation and BH adjustment.
- Deterministic two-gene × 60-sample quickstart.
- Eighteen direct execution tests: 11 compGate tests and 7 joint-model tests.
- All five CMPB figures rebuilt from the 50 packaged Source Data CSVs, with frozen PNG SHA-256 checks.
- Portable benchmark and paired-summary scripts for the full frozen CSVs distributed in the Zenodo candidate.

## Clean installation

Python 3.12 is recommended because the lock file was verified with Python 3.12.10.

```bash
python3.12 -m venv .venv
source .venv/bin/activate          # Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.lock
python -m pip install -e . --no-deps
```

## One-command verification

```bash
python run_checks.py
```

The command runs all 18 tests, executes the quickstart, rebuilds all five figures, verifies their PNG SHA-256 values against the frozen delivery and writes machine-readable results to `verification/local_verification.json`.
It uses only relative paths and the files in this candidate. Matching hashes verify exact PNG reproduction in the verified environment; they are not a claim of byte-identical rendering on every operating system or Matplotlib backend.

## Individual commands

```bash
python tests/test_recap.py
python tests/test_joint_model.py
python scripts/quickstart.py --output verification/quickstart.json
python scripts/rebuild_figure4.py --output verification/figure4
python scripts/rebuild_all_figures.py --output verification/all_figures
```

## Full frozen-data summaries

After downloading the Zenodo candidate, point `--source-data` to its `source_data` directory:

```bash
python analysis/recompute_benchmark_summary.py --source-data /path/to/source_data --output verification/benchmark
python analysis/recompute_paired_summary.py --source-data /path/to/source_data --output verification/paired
```

These two commands recombine frozen per-run/per-cell records; they do not pretend to rerun external cohorts or CIBERSORTx/csSAM.
Full cohort re-estimation requires the public input datasets and schemas in `docs/INPUT_OUTPUT_SCHEMA.md`.

## Scope limits

- The 36-scenario factorial covers seven regression paths; CIBERSORTx and csSAM are separate known-truth arms.
- The 200-replicate real-data analysis is a sample/analysis-unit bootstrap, not donor-level validation.
- Local verification is not independent-team validation. A non-original-developer must execute and sign `docs/REPRODUCTION_RECORD.md`.

See `CHANGELOG.md`, `CITATION.cff`, `docs/INSTALL.md`, `docs/TROUBLESHOOTING.md` and `docs/INPUT_OUTPUT_SCHEMA.md`.
