# Quickstart

The stable simulation entry point is `euskera.evolve`. Model selection keeps
the existing dictionary format:

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
        "gridlength": 4.0, "resol": 32, "tmax": 0.1,
    },
    salva_data_update={"address": "runs/quickstart", "save_number": 2},
)
```

Use a dedicated output directory; generated NPZ files and the `parameters`
file should not be committed. The available model names are `soliton`,
`gaussian_function`, `ell_boson` and `proca`.

Typed settings can be used for evolution, output and diagnostics:

```python
from euskera.evolution import EvolutionConfig
from euskera.io import OutputConfig
from euskera.observables import DiagnosticsConfig

euskera.evolve(
    model_parameters,
    evolution_config=EvolutionConfig(gridlength=4, resol=32, tmax=0.1),
    output_config=OutputConfig(address="runs/quickstart", save_number=2),
    diagnostics_config=DiagnosticsConfig(),
)
```

This does not type model parameters: `model_parameters` remains a dictionary
by design. See [migration](migration.md), [simulations](simulations.md), and
the [notebooks](../examples/).
