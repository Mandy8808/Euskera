# Profiles

## Radial profiles

This directory contains examples for constructing radial profiles that can be used as initial background configurations for simulations.

Open the notebook corresponding to the configuration you want to generate:

- [Linear polarization](example_making_initial_profile_linear.ipynb) — single-field configuration with no self-interaction, $\gamma = 0$.
- [Radial polarization](example_making_initial_profile_radial.ipynb) — single-field configuration with no self-interaction, $\gamma = 1$.
- [Multifrequency](example_making_initial_profile_multifrequency.ipynb) — multifield configuration solved using a shooting method followed by boundary fitting.
- [Self-interaction](example_making_initial_profile_selfinteraction.ipynb) — soliton configuration including a self-interaction term, $\lambda \neq 0$.

For details on the underlying equations and background configurations, see [Background solutions](../../docs/theory/background-solutions.md).

For details on the boundary-fitting procedure used in the multifrequency case, see [Background fitting workflow](../../docs/workflows/background-fitting.md).