# Documentation overview

This directory is the user-facing documentation for Euskera. Read
[`index.md`](index.md) to navigate the guides.

The documentation describes the current public surface without renaming
models or changing numerical behaviour. In particular,
`model_parameters` still uses dictionaries in this release. The typed API
only covers `EvolutionConfig`, `OutputConfig` and `DiagnosticsConfig`.
Initial-model parameters remain model-specific dictionaries, as shown in the
[API guide](api.md) and the [examples](../examples/).

The [architecture](architecture.md), [workflows](workflows.md) and
[scientific background](backgrounds.md) pages provide context; the
[development guide](development.md) explains how to validate documentation
and code changes. The [package map](package-map.md) provides a navigable
overview of the current public modules and symbols.
