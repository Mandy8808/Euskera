# Development and validation

Install the editable package with the test extra:

```bash
python -m pip install -e ".[test]"
```

Run the focused and scientific tests from the repository root:

```bash
python -m pytest
python -m pytest -m scientific
python -m compileall euskera
```

The lightweight Markdown link check is:

```bash
python scripts/check_markdown_links.py
```

It scans Markdown files, ignores external URLs and anchors, and verifies that
relative file targets exist. Documentation changes should preserve the
dictionary-based model API, link to notebooks without moving them, and avoid
committing generated simulation output. The CI workflow is documented in
[`.github/workflows/ci.yml`](../.github/workflows/ci.yml).


## Scientific regression checks

`tests/test_scientific_regressions.py` compares kinetic energy and momentum
against plane-wave formulas, checks the external potential's full weight and
the self-potential's half weight, and compares the assembled linear spectrum
against the union of all three diagonal-sector spectra. Nyquist modes are
included for energy; momentum checks avoid the sign ambiguity at Nyquist.
Eigenvalues are matched without relying on their ordering. Tests use explicit
tolerances and small grids. FFT tests run on NumPy and pyFFTW when installed;
missing pyFFTW is reported as a skip, and the notebook CI job installs it.

`tests/test_run_preparation.py` checks that invalid inputs create no output,
checks physical times and legacy counters in both storage formats, and verifies
metadata and profile copies. These are regression checks, not a substitute
for resolution studies of a physical solution.

## Executable notebook examples

Install the notebook extra and run all three small workflows:

```bash
python -m pip install -e ".[test,notebooks]"
python scripts/run_notebook_smoke.py
# Retain executed notebooks and generated outputs for inspection:
python scripts/run_notebook_smoke.py --output-dir notebook-results
```

The runner executes every cell of each notebook under `examples/workflows/`
in a fresh kernel using the current Python interpreter, with a 120-second
per-cell timeout. `EUSKERA_SMOKE_TEST=1` reduces numerical parameters only.
Each notebook gets a separate temporary working/output directory. The runner
fails on cell errors and never writes execution outputs into source notebooks.
Without `--output-dir`, all generated outputs are temporary. CI runs the same
command and retains executed notebooks and output files as artifacts.

These new self-contained examples complement the existing research notebooks;
the runner does not claim to execute the large research notebooks. The
background example illustrates initial-value profile integration, not a
converged bound-state solution. No video encoder or external datasets are needed.
