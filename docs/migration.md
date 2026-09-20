# Migration and compatibility

This release does **not** migrate the model API. `model_parameters` continues
to be a dictionary mapping existing model names to lists of model-specific
parameter dictionaries:

```python
model_parameters = {"gaussian_function": [{...}]}
```

Keep existing calls using `simulation_parameters_update`,
`salva_data_update` and `comp_conserv_update`; they remain supported.

Typed configuration is intentionally limited to:

- `EvolutionConfig` for numerical evolution settings;
- `OutputConfig` for output location and cadence;
- `DiagnosticsConfig` for conserved quantities and analyses.

These objects may be passed alongside dictionary model parameters. When both
legacy mappings and typed settings are supplied, explicit typed settings win
on overlapping evolution/output/diagnostic keys. No model names or numerical
algorithms changed. See [API](api.md) and [simulations](simulations.md).
