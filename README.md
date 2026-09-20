# Euskera

Euskera is a Python package for numerical evolution of self-gravitating
quantum fields (Schrödinger–Poisson and Gross–Pitaevskii–Poisson), including
soliton, Gaussian, boson-star and Proca-star models.

## Install

Euskera requires Python 3.10 or newer:

```bash
python -m pip install .
# development and tests
python -m pip install -e ".[test]"
```

`pyFFTW` is optional (`python -m pip install ".[fftw]"`); NumPy FFTs are the
fallback. See [installation](docs/installation.md) for details.

## Quickstart

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
    simulation_parameters_update={"gridlength": 4.0, "resol": 32, "tmax": 0.1},
    salva_data_update={"address": "runs/quickstart", "save_number": 2},
)
```

`model_parameters` intentionally remains a dictionary-based API in this
version. Typed configuration is currently limited to `EvolutionConfig`,
`OutputConfig` and `DiagnosticsConfig`; it does not replace model
configuration. See the [quickstart](docs/quickstart.md) and
[migration notes](docs/migration.md).

## Documentation

The [documentation index](docs/index.md) links to installation, API,
architecture, workflows, simulations, visualisation, development and
scientific background. Existing runnable examples remain in
[`examples/`](examples/).

## Project

The project is distributed under the [GNU GPL](LICENSE.txt). Scientific
background and references are in [`references/`](references/).
