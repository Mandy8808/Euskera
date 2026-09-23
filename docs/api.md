# Public API guide

The public API is re-exported from `euskera`. Imports are organized by
workflow and use only the canonical namespaces.

```python
from euskera.evolution import evolve
from euskera.evolution import EvolutionConfig
from euskera.io import OutputConfig
from euskera.observables import DiagnosticsConfig
from euskera.backgrounds import system
from euskera.core import RealGrid
from euskera.observables import massVal
```

## Running a simulation

### `euskera.evolve`

```python
euskera.evolve(
    model_parameters,
    field_components=1,
    salva_data_update=None,
    simulation_parameters_update=None,
    comp_conserv_update=None,
    info=False,
)
```

`model_parameters` maps a model name to a list of update dictionaries. For
example, a Gaussian initial condition can be configured as follows:

```python
import euskera

model_parameters = {
    "gaussian_function": [{
        "positions_gaussiana": [0.0, 0.0, 0.0],
        "amplitude": 1.0,
        "sigma": [1.0, 1.0, 1.0],
    }]
}

euskera.evolve(
    model_parameters,
    simulation_parameters_update={
        "gridlength": 4.0,
        "resol": 64,
        "tmax": 0.5,
    },
    salva_data_update={
        "address": "simulation_data/gaussian_run",
        "save_number": 5,
    },
)
```

For new code, typed configurations validate values at construction time:

```python
evolution = EvolutionConfig(gridlength=4.0, resol=64, tmax=0.5)
output = OutputConfig(address="simulation_data/gaussian_run", save_number=5)
diagnostics = DiagnosticsConfig(Numb_Part=True, Energ=True)
euskera.evolve(model_parameters, evolution_config=evolution,
               output_config=output, diagnostics_config=diagnostics)
```

The canonical configurations take precedence over the compatibility mapping
arguments `simulation_parameters_update`, `salva_data_update` and
`comp_conserv_update` when both are supplied. These mappings remain supported
unchanged. `to_dict()` and `from_legacy()` bridge dataclasses and
the dictionaries consumed by the numerical kernels.

Important simulation options include `lambda_value`, `num_threads`,
`gridlength`, `resol`, `step_factor`, `t0`, `tmax`, `rmax`, `Plim`,
`cmass`, `plott0`, `Boverlap` and `methodEnerg`. The defaults are defined in
`euskera.evolution.evolve`.

## Initial models

`euskera.Models` validates and composes model configurations:

```python
models = euskera.Models(
    gaussian_function=[{
        "positions_gaussiana": [0.0, 0.0, 0.0],
        "amplitude": 1.0,
        "sigma": [1.0, 1.0, 1.0],
    }]
)
grid_data, (psi, rho_i) = models.call_model(
    field_components=1,
    parameters_simulation={"gridlength": 4.0, "resol": 32},
)
```

The lower-level builders `build_soliton`, `build_ell`, `build_proca` and
`build_1d_gaussian` are available for specialised workflows. Use the model
registry when possible so parameter validation remains consistent.

## Grids, potentials and diagnostics

The main numerical helpers are:

- `RealGrid(gridlength=10, resol=128)` for real-space coordinates;
- `KGrid(gridlength=10, resol=128, realspace=False)` for Fourier-space data;
- `Upotential(...)` for the Poisson potential;
- `Npar(psi, Vcell)` for particle number;
- `Energ(...)` and the component energy helpers for energy diagnostics;
- `main_frequency(...)`, `frequMet1(...)` and `frequMet2(...)` for frequency
  analysis.

These functions operate on NumPy arrays produced by the grid and model
stages. Their exact array shapes are best illustrated by the notebooks in
[`examples/`](../examples/).

## Background workflow

`euskera.backgrounds` provides the independent background-solution workflow.
Its
main entry points include:

- `system`, `systemMultifrequency` and `systemMultFreqTot` for background
  systems;
- `shoot` and `freq_shoot` for shooting procedures;
- `profilesFromSolut` for extracting profiles;
- `cheb` and the spectral-operator helpers for Chebyshev discretisation;
- `energEng` and `massVal` for energy and mass calculations.

The canonical modules are `euskera.backgrounds.systems`, `.solvers`,
`.multifrequency`, `.asymptotics`, `.profiles` and `.fitting`:

```python
from euskera.backgrounds import system
system(...)
```

The background API is lower-level than `euskera.evolve`; use the relevant
notebook or source docstring when selecting a solver.

`fitting(..., active=None)` solves the full boundary correction by default.
For identically zero fields, supply a boolean mask in matching unknown and
boundary-equation order, for example `active=[True, True, False]`.
`algebSyst(arg, active=None, info=False)` applies that mask only to the linear
solve; inactive unknowns are zero and the augmented state retains its shape.
See the [fitting workflow](workflows/background-fitting.md) and
[migration notes](migration.md) for replacing `remNul`.

Spectral methods are available from `euskera.spectral`; spectral diagnostic
plots are provided by `euskera.visualization.spectral_plot`.

### Background plotting

`euskera.visualization` provides `plotUsingPerf`, `plotUsingDiscSol`,
`plotUsingSol` and `plotPerf`. They consume discrete profile arrays or
`solve_ivp` results and return the created figure and axes (or axes for plotUsingPerf). Plot settings are
centralized in `euskera.visualization.configuration`.

### Reading HDF5 diagnostics

`euskera.io.read_hdf5_diagnostics(path)` reads schema-1 diagnostic groups and
returns `{snapshot_index: {quantity_name: value}}`, sorted by numeric index.
`StoreSolution(..., diagnostic_names=...)` and
`data_Objgenerator(..., comp_conserv=...)` support named diagnostic writing.
`evolve` supplies the metadata automatically. See
[HDF5 diagnostics](simulations.md#hdf5-diagnostics) for the layout and examples.


### Run validation and timing

`evolve` validates merged configuration and initial fields before creating
output. `t0` must be zero; `tmax` specifies duration. `OutputConfig` accepts
`copy_profiles: bool = False`. `StoreSolution.save_file(data, ti,
physical_time=None)` and `fdata_save(..., physical_time=None)` accept physical
timestamps separately from counters. See [simulations](simulations.md).
