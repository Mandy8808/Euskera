# Changelog

## Unreleased

- Read Gaussian parameters by name and reject incomplete configurations.
- Store HDF5 diagnostics as named numeric datasets, preserving complex samples
  and enabled-quantity metadata; added `read_hdf5_diagnostics` and round-trip tests.

- Preserved sparse Fourier-grid axes in conserved quantities so kinetic energy
  method 2 differentiates along x, y and z correctly. Added analytic plane-wave
  checks for both kinetic methods and multiple field components.

- Replaced right-hand-side zero detection in background fitting with an
  explicit `active` field mask. Updated the background notebook pipeline and
  migration guide; added coupled-system and full augmented-ODE regression tests.

## 1.1.0

- Added reproducible, lightweight smoke coverage for initial profiles, short
  Gaussian workflows, saved-data reading and basic visualisation. Heavy
  notebooks are intentionally not executed in CI.
- Added `scientific` pytest marker and small resolution-comparison tests.
  These checks validate finite results and a minimum resolution consistency;
  they are not a substitute for a production convergence study.
- Promoted the typed `EvolutionConfig`, `OutputConfig` and
  `DiagnosticsConfig` workflow API and documented the canonical imports.
- Legacy update dictionaries remain supported for migration, but new code
  should use typed configurations. No legacy directories or scientific data
  were removed.
