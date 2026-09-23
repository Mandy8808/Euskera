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

## Explicit active fields in background fitting

`fitting` now accepts `active=None` as its final optional argument and passes
it to `algebSyst`. The latter replaces `remNul` with `active`: replace
`remNul=False` with `active=None`; replace `remNul=True` with a boolean mask
such as `active=[True, True, False]`, based on the physical configuration.
Update positional calls as well; a single boolean is not an active-field mask.

Omitting the mask now solves the complete correction system, even when the
right-hand side contains zeros. Configurations with inactive fields must pass
their mask explicitly to avoid a singular full system. The augmented ODE state
and downstream profile format are unchanged. See the
[fitting workflow](workflows/background-fitting.md) for the notebook pipeline.

## Gaussian inputs and HDF5 diagnostics

Gaussian dictionary keys may appear in any order. Every Gaussian must explicitly
supply position, amplitude and widths; incomplete entries now raise `ValueError`.

HDF5 diagnostics now use named numeric datasets inside versioned snapshot groups.
Use `read_hdf5_diagnostics` instead of treating `save_energies_<index>` as a single
array. Direct HDF5 diagnostic writers must supply `diagnostic_names`; the
simulation pipeline does this automatically. NPZ diagnostics are unchanged.
See [simulation output](simulations.md#hdf5-diagnostics).


## Validated runs and physical time

Nonzero `t0` is now rejected rather than silently ignored. `tmax` remains the
simulation duration starting at zero. Unknown saving flags/formats, nonfinite
parameters and invalid model configurations fail before output is created.
`OutputConfig(copy_profiles=True)` optionally archives input radial arrays.

Saved `t` retains its legacy counter meaning. Prefer the new `time` dataset
for physical times and `snapshot_index` for indices. Direct storage callers
must pass `physical_time` to obtain a `time` dataset. Existing NPZ/HDF5 readers
can continue using their old keys; readers that enumerate every archive entry
as a snapshot should filter by snapshot name and exclude metadata datasets.
Every new simulation also writes `run_metadata.json`; see
[run provenance](simulations.md#run-provenance-and-optional-profile-copies).
