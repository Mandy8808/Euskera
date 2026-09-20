<p align="center">
  <img src="galleries/encabezado.png" width="760" alt="Euskera">
</p>

<p align="center">
  <a href="LICENSE.txt"><img src="https://img.shields.io/badge/license-GNU%20GPL-green" alt="GNU GPL license"></a>
  <a href="https://www.python.org"><img src="https://img.shields.io/badge/python-3.10%2B-blue" alt="Python 3.10 or newer"></a>
  <a href="https://github.com/Mandy8808/Euskera/actions/workflows/ci.yml"><img src="https://github.com/Mandy8808/Euskera/actions/workflows/ci.yml/badge.svg" alt="Continuous integration"></a>
</p>

# Euskera

Euskera is a Python code for constructing and evolving self-gravitating
quantum fields. It contains numerical tools for the
[Schrödinger-Poisson (SP)](#theoretical-scope) and
Gross-Pitaevskii-Poisson (GPP) systems, together with models for solitons,
boson stars and Proca stars.

The repository is both a Python distribution and a collection of scientific
resources. Use the links below to navigate directly to the part you need.

## Navigation

- [Theoretical scope](#theoretical-scope)
- [Repository map](#repository-map)
- [Installation](#installation)
- [First steps](#first-steps)
- [Examples and scientific resources](#examples-and-scientific-resources)
- [Technical documentation](#technical-documentation)
- [Development and validation](#development-and-validation)
- [Citation and contact](#citation-and-contact)

## Theoretical scope

The SP system describes a non-relativistic, self-gravitating scalar field
without self-interaction. In schematic form,

```text
i dPsi/dt = -(1/2) Laplacian(Psi) + U Psi
Laplacian(U) = 4 pi G rho
```

Euskera extends this framework with a nonlinear self-interaction term, giving
the GPP system. The same numerical evolution machinery can be initialized
with different field configurations, including:

- soliton profiles;
- Gaussian perturbations;
- ellipsoidal boson-star configurations;
- Proca-star configurations, including multi-frequency cases.

The evolution uses a split-step Fourier method. Spatial grids, potentials,
conserved quantities, profile generation and saved-output visualisation are
exposed through the package API. The mathematical background and numerical
derivation are collected in [`references/`](references/), especially
[`references/main.tex`](references/main.tex) and
[`references/SP_system.pdf`](references/SP_system.pdf).

## Repository map

The distribution is the union of two importable Python packages:

| Path | Role |
| --- | --- |
| [`euskera/`](euskera/) | Main evolution code, field models, grids, potentials, plotting and data saving |
| [`background/`](background/) | Background solutions, spectral methods and radial-profile tools |

The rest of the repository contains supporting scientific material:

| Path | Contents |
| --- | --- |
| [`examples/`](examples/) | Jupyter notebooks for initialization, evolution and visualisation |
| [`profiles/`](profiles/) | Numerical radial profiles used by the examples |
| [`simulation_data/`](simulation_data/) | Representative saved simulation outputs |
| [`galleries/`](galleries/) | Figures and example videos |
| [`references/`](references/) | Equations, bibliography and source material |
| [`tests/`](tests/) | Import, initialization and fallback tests |
| [`pyproject.toml`](pyproject.toml) | Package metadata and dependency declarations |
| [`.github/workflows/ci.yml`](.github/workflows/ci.yml) | Continuous-integration workflow |

## Installation

Euskera supports Python 3.10 and newer.

Install the package from a local clone:

```bash
python -m pip install .
```

For development and tests:

```bash
python -m pip install -e ".[test]"
```

The optional [pyFFTW](https://pyfftw.readthedocs.io/) acceleration can be
installed with:

```bash
python -m pip install ".[fftw]"
```

`pyFFTW` is not required. Without it, Fourier transforms use NumPy. FFmpeg is
also optional and is only needed when exporting animations to video.

The package has two public top-level imports:

```python
import euskera
from euskera.backgrounds import system
```

## First steps

The main evolution entry point is `euskera.evolve`. The model registry
currently provides:

```text
soliton
gaussian_function
ell_boson
proca
```

The notebooks in [`examples/`](examples/) are the best source for complete
parameter configurations. A minimal API smoke check is:

```python
import euskera
from euskera.backgrounds import system

print(euskera.Models)
print(euskera.evolve)
print(system)
```

Simulation output is written as NPZ by default. HDF5 output is available
through the saving API when `h5py` is installed. Generated output should
normally be written to a local directory rather than added to version
control.

## Examples and scientific resources

The examples cover:

- construction of initial profiles;
- Gaussian, soliton, boson-star and Proca-star simulations;
- collisions and multi-frequency configurations;
- extraction and visualisation of saved data;
- frequency and oscillation analysis.

The repository also includes published profiles, representative simulation
data and videos so that results can be inspected without running a complete
simulation. The scientific description is available in
[`references/main.tex`](references/main.tex), with bibliography data in
[`references/bibliografía.bib`](references/bibliografía.bib).

## Technical documentation

The [`docs/`](docs/) directory contains the technical documentation:

- [Documentation index](docs/index.md): orientation and links to all guides;
- [Architecture](docs/architecture.md): package boundaries and the simulation
  data flow;
- [API guide](docs/api.md): public entry points, models and numerical helpers;
- [Workflows](docs/workflows.md): practical recipes for configuring,
  executing and inspecting simulations.

## Development and validation

Run the test suite:

```bash
python -m pytest
```

Compile both packages:

```bash
python -m compileall euskera
```

Build a wheel:

```bash
python -m pip wheel . --no-deps
```

The GitHub Actions workflow tests the package on Python 3.10, 3.11, 3.12,
3.13 and 3.14. It installs the package, compiles the sources and runs the
tests.

## Citation and contact

The project is distributed under the GNU GPL license; see
[`LICENSE.txt`](LICENSE.txt).

If Euskera contributes to a project that leads to a publication, please
acknowledge the authors and cite the associated scientific work.

Contact:

- alberto.diez(at)fisica.ugto.mx
- arestrada(at)fisica.uaz.edu.mx
- arestrada(at)fisica.ugto.mx
