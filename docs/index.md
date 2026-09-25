# Euskera documentation

Start with the [repository README](../README.md) for installation and the
shortest possible example. This index groups the longer guides by task.

## Getting started

- [Installation](installation.md): dependencies, optional acceleration and
  development setup.
- [Quickstart](quickstart.md): the smallest dictionary-based simulation and
  how to choose an example notebook.
- [Migration](migration.md): compatibility guidance and the deliberately
  limited typed configuration API.

## Understanding the package

- [Package map](package-map.md): navigable modules and public symbols.
- [Background overview](backgrounds.md): physical context and background
  module overview.
- [Architecture](architecture.md): package boundaries and simulation data
  flow.
- [API guide](api.md): public functions, models, grids and diagnostics.
- [Simulations](simulations.md): reproducible run configuration and output.
- [Visualisation](visualization.md): inspecting saved data and exporting
  figures or videos.
- [Practical workflows](workflows.md): notebook-oriented usage.

## Theory and methods

- [Background solutions](theory/background-solutions.md): equations and
  radial background families from arXiv:2412.06901.
- [Spectral method](theory/spectral-method.md): Chebyshev discretisation and
  stability analysis from arXiv:2512.04376.
- [Spectral stability notebook](../examples/spectral/example_SprectralStability_proca.ipynb):
  executable Proca background spectrum example.
- [Boundary-value method](methods/boundary-value-method.md): numerical
  parameter correction from arXiv:2208.13221.

## Workflows

- [Background shooting](workflows/background-shooting.md): node and frequency
  selection for radial solutions.
- [Background fitting](workflows/background-fitting.md): independent
  boundary-condition parameter correction.
- [Dynamical simulation](workflows/simulation.md): three-dimensional time
  evolution and diagnostics.

## Contributing

- [Development](development.md): tests, link checking and contribution
  conventions.

## Repository resources

- [Examples and notebooks](../examples/)
- [Profiles](../profiles/)
- [Representative simulation data](../simulation_data/)
- [Scientific references](../references/)

- [Executable workflow notebooks](../examples/workflows/README.md)
