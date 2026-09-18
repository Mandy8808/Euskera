<p align="center" width="160%">
    <img width="80%" src="galleries/encabezado.png" alt="Euskera">
</p>

<p float="left">
<a href="LICENSE.txt"><img src="https://img.shields.io/badge/GNU%20GPL%20license-green" alt="GNU GPL license"></a>
<a href="https://www.python.org"><img src="https://img.shields.io/badge/Language-Python-blue" alt="Python"></a>
<a href="https://github.com/Mandy8808/Euskera/actions/workflows/ci.yml"><img src="https://github.com/Mandy8808/Euskera/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
</p>

# Euskera

Euskera is a Python package for the numerical time evolution of the
Schrodinger-Poisson and Gross-Pitaevskii-Poisson systems. It is based on the
public code [PyUltraLight](https://github.com/auckland-cosmo/PyUltraLight).

The project extends the original approach by adding a self-interaction term
and support for multi-frequency Proca stars and boson-star models. The
mathematical and numerical procedure is described in the reference material
included in [`references/`](references/).

The distribution combines two related packages:

- [`euskera/`](euskera/): time evolution, physical models, grids, potentials,
  conserved quantities, plotting and data saving.
- [`background/`](background/): background solutions, spectral methods,
  profile construction and related plotting utilities.

## Project contents

| Path | Description |
| --- | --- |
| [`euskera/`](euskera/) | Main simulation and modelling package |
| [`background/`](background/) | Background and spectral-method package |
| [`examples/`](examples/) | Jupyter notebooks showing representative workflows |
| [`profiles/`](profiles/) | Published numerical radial profiles |
| [`simulation_data/`](simulation_data/) | Example simulation outputs |
| [`galleries/`](galleries/) | Figures and example videos |
| [`references/`](references/) | Scientific reference material and bibliography |
| [`tests/`](tests/) | Automated package and smoke tests |

## Requirements

Euskera requires Python 3.10 or newer.

The core dependencies are:

- [NumPy](https://numpy.org)
- [SciPy](https://scipy.org)
- [NumExpr](https://numexpr.readthedocs.io/)
- [Numba](https://numba.pydata.org/)
- [Matplotlib](https://matplotlib.org/)
- [Cycler](https://cycler.readthedocs.io/)
- [h5py](https://www.h5py.org/)
- [pandas](https://pandas.pydata.org/)
- [IPython](https://ipython.org/)

[`pyFFTW`](https://pyfftw.readthedocs.io/) is optional. When installed, it
provides accelerated Fourier transforms. If it is not installed, Euskera uses
the NumPy FFT implementation instead.

[FFmpeg](https://ffmpeg.org/) is also optional and is only needed to export
animations to video.

## Installation

### From a local clone

```bash
python -m pip install .
```

For development and testing:

```bash
python -m pip install ".[test]"
```

To enable the optional FFTW acceleration:

```bash
python -m pip install ".[fftw]"
```

After installation, both parts of the distribution are importable:

```python
import euskera
import background
```

### Development environment

The repository contains the package metadata in
[`pyproject.toml`](pyproject.toml). An editable installation is useful when
developing:

```bash
python -m pip install -e ".[test]"
```

The repository CI tests Python 3.10 through 3.14 and verifies installation,
compilation and the automated tests.

## Basic usage

The main entry point for a time evolution is `euskera.evolve`. Models are
selected by name and configured through parameter dictionaries. For example,
the public API can be inspected with:

```python
import euskera

print(euskera.Models)
print(euskera.evolve)
```

Available model names are:

- `soliton`
- `gaussian_function`
- `ell_boson`
- `proca`

The notebooks in [`examples/`](examples/) contain complete configurations for
the supported models, profile generation, collision simulations and
visualisation of saved data. They should be used as the authoritative examples
for the detailed parameter shapes.

By default, simulation output is saved in NPZ format. HDF5 output is available
through the data-saving API when `h5py` is installed. Simulation output should
normally be written to a local output directory rather than committed to the
repository.

## Validation and development commands

Run the tests:

```bash
python -m pytest
```

Compile the Python sources:

```bash
python -m compileall euskera background
```

Build a wheel:

```bash
python -m pip wheel . --no-deps
```

## Scientific material and citation

The [`references/`](references/) directory contains the current reference
material, including the system description, bibliography and source document.
The repository also includes published profiles and representative simulation
outputs for comparison and visualisation.

If Euskera contributes to a project that leads to a publication, please
acknowledge the authors and cite the associated work.

## License

This project is distributed under the GNU GPL license. See
[`LICENSE.txt`](LICENSE.txt).

## Contact

- alberto.diez(at)fisica.ugto.mx
- arestrada(at)fisica.uaz.edu.mx
- arestrada(at)fisica.ugto.mx
