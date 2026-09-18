<p align="center" width="160%">
    <img width="80%" src="galleries/encabezado.png">
</p>

<p float="left">
<a href="LICENSE.txt"><img src="https://img.shields.io/badge/GNU%20GPL%20license-green" alt="GNU GPL license"></a>
<a href="https://www.python.org"><img src="https://img.shields.io/badge/Language-Python-blue" alt="Python"></a>
</p>

# Euskera

Euskera is a Python package for the numerical time evolution of
Schrodinger-Poisson and Gross-Pitaevskii-Poisson systems. It extends the
PyUltraLight approach with self-interaction, multi-frequency Proca stars and
boson-star models.

The distribution contains two related namespaces:

- [`euskera/`](euskera/): time evolution, physical models, plotting and data
  saving.
- [`background/`](background/): background and spectral methods used to build
  numerical profiles.

## Requirements

- Python 3.10 or newer.
- NumPy, SciPy, NumExpr and Numba for the numerical core.
- Matplotlib and Cycler for plotting.
- h5py for HDF5 output.
- pandas for spectral-method utilities.
- IPython for video visualisation helpers.

`pyFFTW` is optional. When it is installed, Euskera uses it for accelerated
Fourier transforms; otherwise it falls back to NumPy's FFT implementation.

FFmpeg is optional and is only needed when exporting generated animations to
video.

## Installation

Install the package from a clone of this repository:

```bash
python -m pip install .
```

For development and tests:

```bash
python -m pip install ".[test]"
```

To enable the optional FFTW acceleration:

```bash
python -m pip install ".[fftw]"
```

The package can then be imported normally:

```python
import euskera
import background
```

## Quick checks

Run the test suite with:

```bash
python -m pytest
```

Compile all package sources with:

```bash
python -m compileall euskera background
```

The repository includes example notebooks in [`examples/`](examples/),
published profiles in [`profiles/`](profiles/), simulation data in
[`simulation_data/`](simulation_data/) and example videos in
[`galleries/Videos/`](galleries/Videos/).

By default, simulation output is saved in NPZ format. HDF5 output is also
available through the data-saving API when h5py is installed.

## Citation and contact

If Euskera contributes to a project that leads to a publication, please
acknowledge the authors and cite the associated work.

Contact:

- alberto.diez(at)fisica.ugto.mx
- arestrada(at)fisica.uaz.edu.mx
- arestrada(at)fisica.ugto.mx
