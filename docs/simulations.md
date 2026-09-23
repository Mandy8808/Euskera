# Running simulations

Configure the physical model with dictionary-based `model_parameters`, then
choose numerical evolution and output settings. The four model names are
documented in [architecture](architecture.md) and their parameter shapes are
shown in the [API guide](api.md) and [examples](../examples/).

For new runs, explicit typed settings make the three supported configuration
groups clear:

```python
from euskera.evolution import EvolutionConfig
from euskera.io import OutputConfig
from euskera.observables import DiagnosticsConfig

evolution = EvolutionConfig(gridlength=10, resol=128, tmax=1)
output = OutputConfig(address="runs/my_experiment", save_number=20)
diagnostics = DiagnosticsConfig()
euskera.evolve(model_parameters, evolution_config=evolution,
               output_config=output, diagnostics_config=diagnostics)
```

The legacy mappings `simulation_parameters_update`, `salva_data_update` and
`comp_conserv_update` remain supported. NPZ snapshots and a parameters file
are written under `address` by default; HDF5 is available through the saving
API when `h5py` is installed. Keep run outputs outside tracked source
directories and record the parameter file with results.

For reproducibility, keep the notebook, configuration and package version
together. See [workflows](workflows.md) and [visualisation](visualization.md).

## Gaussian configurations

Each entry in `gaussian_function` must provide `positions_gaussiana`,
`amplitude`, and `sigma`. Dictionary key order has no effect:

```python
model_parameters = {
    "gaussian_function": [{
        "amplitude": 1.0,
        "sigma": [1.0, 0.8, 1.2],
        "positions_gaussiana": [0.0, 0.0, 0.0],
    }]
}
```

Missing parameters raise `ValueError` instead of silently dropping a Gaussian.
Direct users of `Gaussiana_Model` must supply equally long parameter lists.

## HDF5 diagnostics

Select `OutputConfig(format="hdf5", ...)` to write HDF5 output. When
`save_energies` is enabled, `evolve` passes the enabled diagnostic names to
storage automatically. The consolidated `end_save_energies.hdf5` contains
one group per saved snapshot:

```text
save_energies_0/
    mass_total
    mass_components
    energy
    momentum
    frequency_samples
save_energies_1/
    ...
t
```

Only enabled quantities are present. `Numb_Part` produces `mass_total` and
`mass_components`; `Energ` produces `energy`; `Pi` produces the three-vector
`momentum`. `Frequency` produces complex `frequency_samples`: these are
wavefunction values at the tracked positions, not fitted frequencies.
Each snapshot group records `schema_version=1` and `diagnostic_names` in
calculation order. With all diagnostics disabled, snapshot groups are empty.

The group suffixes and legacy `t` dataset remain save counters. New runs also
save `snapshot_index` (explicit counters) and `time` (physical times in code
units), supplied by the integrator at each completed save step. Use `time`
for analysis instead of reconstructing it from the number of files. Current
runs require `t0=0`; `tmax` is duration, and nonzero `t0` is rejected.
Direct `StoreSolution.save_file` calls may pass `physical_time`; the `time`
dataset is omitted when any timestamp was not supplied, rather than guessed.

Read the named quantities with:

```python
from euskera.io import read_hdf5_diagnostics

snapshots = read_hdf5_diagnostics("runs/my_experiment/end_save_energies.hdf5")
for counter, quantities in snapshots.items():
    print(counter, quantities.get("energy"))
```

The reader returns a dictionary in numerical snapshot order; each value is a
quantity dictionary containing numeric scalars or arrays. It checks the schema
version and does not read legacy unlabelled diagnostic datasets. Standard HDF5
field and grid datasets remain readable with `h5py` as before. NPZ output and
`Conserv`'s positional return format are unchanged.

For direct storage calls, provide names in the exact order of the `Conserv`
output, for example:

```python
from euskera.io import StoreSolution

store = StoreSolution(
    "runs/my_experiment", "save_energies", format="hdf5",
    diagnostic_names=["Numb_Part", "Energ"],
)
# cData must contain [particle_number_result, energy_result] in that order.
store.save_file(cData, ti=0)
store.close_file("end_save_energies")
```

Alternatively, pass `comp_conserv` to `data_Objgenerator`. Diagnostic names are
required for HDF5 diagnostic storage; the writer does not guess them from
array positions or shapes.


## Input validation and run preparation

Both legacy dictionaries and typed configurations are validated after merging,
before output directories or parameter files are created. Initial fields,
potential and enabled initial energy are constructed and checked before output.
Numeric settings must be finite; Gaussian widths, profile spacing and scales
must be positive. Each radial profile must be a finite real one-dimensional
array with at least two samples. Configurations require consistent component
counts, three components for Proca polarization, and `2*ell+1` for ell models.
Unknown output formats, saving flags and diagnostic options are rejected.
Only energy methods 1 and 2 are accepted. A supplied dataclass that was mutated
after construction is revalidated at run start.

## Run provenance and optional profile copies

Every `evolve` run writes `run_metadata.json` alongside `parameters.txt`. It
records a unique run ID, creation time, Python and package versions, FFT backend,
effective model/numerical/output/diagnostic settings, integration step size and
step count. Effective model scales reflect the existing alpha=1 rule when
self-interaction is enabled. Git revision and dirty state are recorded when
available; unavailable values are null, not evidence of a clean checkout.
A dirty flag does not archive uncommitted source changes.

Radial profile descriptors contain SHA-256 hashes computed from dtype, shape
and contiguous array bytes. Hashes identify data but cannot recreate them.
To save the actual arrays as numeric entries in `input_profiles.npz`, use:

```python
output = OutputConfig(address="runs/reproducible", copy_profiles=True)
```

The profile descriptor's `copy_key` identifies its entry in that archive.
Gaussian-only runs have no input radial profiles to copy. Keep the metadata,
profile archive and matching source revision together for reproducibility.
NPZ and HDF5 snapshots both retain `t` and add `snapshot_index` and `time`.
