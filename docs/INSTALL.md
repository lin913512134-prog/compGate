# Installation and verification

1. Use a clean Python 3.12 environment on Linux or the environment actually used by the authors.
2. Install `requirements.lock`, then install this package editable with `--no-deps`.
3. Run `python run_checks.py` from the repository root.
4. Confirm 18/18 tests, zero quickstart failures and successful Figure 4 reconstruction.
5. Record OS, Python, package version, commit/tag, runtime, memory, warnings and failures in `docs/REPRODUCTION_RECORD.md`.

No step requires the original author's local directory structure. Do not use a different package version and then report it as this candidate's reproduction.
