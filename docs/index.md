# Euskera technical documentation

This section describes the package from the perspective of a user who needs
to configure simulations, extend an initial model or inspect saved results.
For a quick introduction, start with the [repository README](../README.md).

## Guides

- [Architecture](architecture.md) explains how grids, models, potentials,
  evolution and output saving interact.
- [API guide](api.md) lists the public entry points exported by `euskera` and
  its workflow namespaces.
- [Workflows](workflows.md) contains practical configuration examples and
  guidance for notebooks, output and optional acceleration.

## Scope and conventions

Euskera evolves discretised self-gravitating fields on a three-dimensional
Cartesian grid. Unless stated otherwise:

- `gridlength` is the physical length of the cubic domain;
- `resol` is the number of points per spatial direction;
- arrays use NumPy conventions and may contain one leading axis for field
  components;
- simulation configuration is passed as dictionaries so that notebooks can
  be adapted without changing package code.

The notebooks in [`examples/`](../examples/) remain the most complete
reference for scientific parameter combinations. These pages document the
stable public surface, not every internal implementation detail.
