# Installation

Euskera supports Python 3.10 or newer. From a clone of the repository:

```bash
python -m pip install .
```

For an editable development installation with tests:

```bash
python -m pip install -e ".[test]"
```

The runtime dependencies are declared in [`pyproject.toml`](../pyproject.toml).
NumPy FFTs are always available. Install the optional `pyFFTW` extra for FFT
acceleration:

```bash
python -m pip install ".[fftw]"
```

FFmpeg is only needed when exporting animations. HDF5 output uses the
optional `h5py` dependency already declared by the package.

Verify the installation with:

```bash
python -c "import euskera; print(euskera.evolve)"
```

Continue with the [quickstart](quickstart.md), or install the test extra and
read [development](development.md).
