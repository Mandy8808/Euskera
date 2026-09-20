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
