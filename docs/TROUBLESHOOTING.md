# Troubleshooting

- **`No module named compgate`:** activate the clean environment and run `python -m pip install -e . --no-deps` from the repository root.
- **Dependency conflict:** recreate the environment and install `requirements.lock` before installing the package.
- **Matplotlib font differences:** use the locked matplotlib version and record the OS/font environment; numerical anchors are authoritative even if font rasterization differs.
- **Figure inputs missing:** run from the repository root and preserve `examples/figure4_source/`; do not rename the five CSVs.
- **Full frozen-data inputs missing:** download the Zenodo `source_data` directory and pass its path explicitly to the analysis scripts.
- **Do not interpret scope literally beyond the evidence:** seven regression paths are not every SoA method; analysis-unit bootstrap is not donor-level validation.
