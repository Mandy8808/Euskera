# Changelog

## Unreleased

- Added analytic spectral, kinetic, momentum and central-potential regression checks.
- Validate merged configurations and initial fields before creating run output.
- Record physical snapshot times, run metadata and optional hashed profile copies.
- Added three self-contained workflow notebooks, a fresh-kernel smoke runner and CI.


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

## Selective output (unreleased)

- Temporal rules by variable and geometry; three planes and axes for rho, psi, and phi.
- Independent diagnostic calendars, potential half-step synchronization, and preserved global indices.
- NPZ/HDF5 reading and plots with coordinates, orientation, and physical time; data_save compatibility.
- Selective output guide and notebooks, numerical tests, and a change report.

## Consolidation integration (unreleased)

- Integrated recoverable NPZ/HDF5 consolidation with after-success and incremental cleanup.
- Applied cleanup options to legacy and selective output, preserving original sample times during recovery.
- Fixed serialized output rules in legacy configuration, broken profile documentation links, and CI workflow-notebook LFS downloads.
