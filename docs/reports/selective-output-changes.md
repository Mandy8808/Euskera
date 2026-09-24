# Implementation report: selective output

Historical report: see [repository reconciliation](reconciliation.md) for the subsequently integrated recoverable consolidation.

Date: 2026-09-24. Documentation and newly added notebook content have been translated into English.

## Result

- Independent rules by variable (`rho`, `psi`, `phi`) and geometry: volumes, xy/xz/yz planes, and x/y/z lines.
- `All`, `Last`, `TimeRange`, and `Final` selections, subsampling with `every`, and optional initial-state inclusion per rule.
- A common calendar with global indices and physical times; overlapping rules do not duplicate samples.
- Independent diagnostics: particle number, energy, momentum, and frequency samples. Their requested union is computed once per instant.
- The pending potential half-step is completed once per requested instant and at termination when needed. FFTs and potential updates continue at every integration step.
- NPZ and HDF5, spatial metadata, a common reader, slice plots, and temporal alignment in the existing animation.
- Compatibility: `rules=None` retains `data_save`; explicit rules replace the existing output selection.

## Verification

- Complete test suite in **ciencia**, with **pyFFTW 0.15.1**: **144 passed, none skipped**.
- NumPy was also checked in `.venv`: the initial suite had 135 passes and 7 skips due to unavailable FFTW, followed by 8 passes in the expanded selective-output module.
- Numerical comparisons: partial versus full output for rho, complex psi, phi, and all four diagnostics; also complete diagnostics with partial field output.
- Spatial tests: all three planes and lines for each variable, two psi components, and coordinates without a node exactly at zero.
- Both formats: NPZ and HDF5; empty windows, overlapping rules, serialized configuration, and validation errors.
- Four workflow notebooks executed successfully: background_profile, gaussian_simulation, selective_output, and visualize_generated_data. Verification artifacts are in `/tmp/euskera-selective-notebook-check` (temporary storage).
- Eleven new simulation-notebook configuration cells were checked for valid syntax. Their long simulations were not executed and their historical results were not regenerated.
- The package map and Markdown links were checked. The subsequent English-only edits do not change numerical behavior; notebook syntax and documentation links were checked again.

## Documentation and notebooks

The [selective output guide](../selective-output.md) defines parameters, slices, calendars, half-steps, diagnostics, reading, and limitations. The [new notebook](../../examples/workflows/selective_output.ipynb) runs and verifies a small example.

Simulation notebooks include optional configuration before each evolution call: `USE_SELECTIVE_OUTPUT=False` preserves their existing analysis; setting it to `True` enables rules and a directory with the `_selective` suffix. A new directory is recommended for each run. These examples use `Last(min(100, save_number + 1))` to avoid requesting more samples than are available; a final window of 100 samples with earlier history requires a larger calendar.

## Compatibility and deliberate limitations

- `rules=None` keeps the existing writers, names, initial sample, and common output cadence. An explicit list uses the new scheduler instead; `rules=[]` writes only grid and metadata.
- `output_config` takes precedence over `salva_data_update` when both are provided, as before. New output names require the new reader or adapted analysis code.
- The existing time-step adjustment involving `save_number` is unchanged.
- Adaptive termination and rolling retention for unknown endpoints are not implemented.
- Slices use the node nearest to zero, without interpolation or arbitrary positions.
- Initial validation diagnostics are retained even when frame 0 is not saved.
- `read_output` loads a complete output into memory; direct per-sample HDF5 access is recommended for large volumes.
- NPZ diagnostics retain object encoding and require trusted files; HDF5 uses numeric datasets.
- Small rounding differences from combining half-steps are checked with numerical tolerances.

## Line-report scope

The tables compare against a snapshot taken **at the start of this implementation**, not against HEAD. They exclude pre-existing user changes in potential, observables, documentation, and notebooks. Line numbers are one-based file positions; notebook numbers refer to JSON lines, not cell numbers. The [detailed patch](selective-output.patch) records exact additions and deletions, including new files. This report and the patch are excluded from their own comparison. Existing source text shown as deleted patch lines is reproduced verbatim.

### `euskera/io/schedule.py`

+108 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 0 | 1–108 |

### `euskera/io/selected_output.py`

+123 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 0 | 1–123 |

### `euskera/io/config.py`

+10 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 4 | 5–5 |
| insertion after 23 | 25–25 |
| insertion after 25 | 28–35 |

### `euskera/io/save_data.py`

+1 / −1 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| 224–224 | 224–224 |

### `euskera/io/__init__.py`

+4 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 8 | 9–12 |

### `euskera/__init__.py`

+4 / −1 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| 85–85 | 85–88 |

### `euskera/evolution/evolve.py`

+14 / −3 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 25 | 26–26 |
| insertion after 145 | 147–149 |
| 206–206 | 210–211 |
| 210–210 | 215–221 |
| 224–224 | 235–235 |

### `euskera/evolution/evolut_routines.py`

+41 / −18 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| 25–25 | 25–25 |
| 82–88 | 82–108 |
| 90–96 | deletion after 109 |
| 107–108 | 120–127 |
| 112–112 | 131–135 |

### `euskera/visualization/selected_output.py`

+37 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 0 | 1–37 |

### `euskera/visualization/__init__.py`

+1 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 12 | 13–13 |

### `euskera/visualization/video_make.py`

+27 / −2 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 239 | 240–264 |
| 242–242 | 267–267 |
| 248–248 | 273–273 |

### `tests/test_output_rules.py`

+121 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 0 | 1–121 |

### `docs/selective-output.md`

+110 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 0 | 1–110 |

### `docs/simulations.md`

+4 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 154 | 155–158 |

### `docs/api.md`

+4 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 176 | 177–180 |

### `docs/visualization.md`

+4 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 27 | 28–31 |

### `docs/workflows/simulation.md`

+4 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 48 | 49–52 |

### `docs/package-map.md`

+17 / −7 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| 54–54 | 54–54 |
| 60–60 | 60–60 |
| insertion after 73 | 74–81 |
| 122–126 | 130–134 |
| insertion after 180 | 189–190 |

### `README.md`

+2 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 57 | 58–59 |

### `CHANGELOG.md`

+7 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 35 | 36–42 |

### `examples/workflows/selective_output.ipynb`

+119 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 0 | 1–119 |

### `examples/simulations/example_computing_oscilation_soliton_frequency.ipynb`

+10 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 569 | 570–579 |

### `examples/simulations/example_simulation_colision_proca_radial_linear.ipynb`

+39 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 463 | 464–473 |
| insertion after 467 | 478–505 |
| insertion after 526 | 565–565 |

### `examples/simulations/example_simulation_ell_boson.ipynb`

+39 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 395 | 396–405 |
| insertion after 399 | 410–437 |
| insertion after 458 | 497–497 |

### `examples/simulations/example_simulation_gaussiana.ipynb`

+39 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 180 | 181–190 |
| insertion after 181 | 192–219 |
| insertion after 239 | 278–278 |

### `examples/simulations/example_simulation_multifrequency.ipynb`

+39 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 405 | 406–415 |
| insertion after 409 | 420–447 |
| insertion after 434 | 473–473 |

### `examples/simulations/example_simulation_proca_circular.ipynb`

+39 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 381 | 382–391 |
| insertion after 385 | 396–423 |
| insertion after 410 | 449–449 |

### `examples/simulations/example_simulation_proca_linear.ipynb`

+39 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 397 | 398–407 |
| insertion after 401 | 412–439 |
| insertion after 426 | 465–465 |

### `examples/simulations/example_simulation_proca_radial.ipynb`

+39 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 381 | 382–391 |
| insertion after 385 | 396–423 |
| insertion after 410 | 449–449 |

### `examples/simulations/example_simulation_soliton_gaussiana.ipynb`

+39 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 379 | 380–389 |
| insertion after 383 | 394–421 |
| insertion after 441 | 480–480 |

### `examples/simulations/example_simulation_soliton_halo.ipynb`

+39 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 397 | 398–407 |
| insertion after 401 | 412–439 |
| insertion after 469 | 508–508 |

### `examples/simulations/example_simulation_soliton_lam_0.ipynb`

+78 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 395 | 396–405 |
| insertion after 399 | 410–432 |
| insertion after 400 | 434–438 |
| insertion after 408 | 447–447 |
| insertion after 644 | 684–693 |
| insertion after 645 | 695–722 |
| insertion after 704 | 782–782 |

### `examples/simulations/example_visualizing_simulation_data.ipynb`

+10 / −0 lines relative to the start of this implementation.

| Previous lines | Current lines |
| --- | --- |
| insertion after 1351222 | 1351223–1351232 |

