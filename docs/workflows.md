# Practical workflows

## Start from a notebook

The recommended workflow is to copy one of the notebooks in
[`examples/`](../examples/) and change only the physical and numerical
parameters needed for the experiment. This preserves the expected array
shapes and output conventions.

The notebooks are grouped around:

- profile construction;
- Gaussian, soliton, ellipsoidal and Proca simulations;
- collisions and multi-frequency configurations;
- saved-data visualisation;
- frequency and oscillation analysis.

## Configure output explicitly

For reproducible runs, set the output address and save cadence instead of
relying on the defaults:

```python
save_config = {
    "format": "npz",
    "address": "runs/my_experiment",
    "save_number": 20,
    "data_save": {
        "grid": True,
        "save_rho": True,
        "save_psi": False,
        "save_phi": False,
        "save_energies": True,
    },
}
```

Pass this mapping as `salva_data_update` to `euskera.evolve`. Keep the
configuration and generated data together so a run can be reconstructed.

## Select conserved quantities

The default diagnostics include particle number and energy. Additional
diagnostics can be enabled with `comp_conserv_update`:

```python
diagnostics = {
    "Numb_Part": True,
    "Energ": True,
    "Pi": True,
    "Ji": False,
    "Frequency": True,
}
```

Only enable expensive diagnostics when they are needed for the analysis,
especially at high resolution.

## Optional Fourier acceleration

NumPy FFTs are the default fallback. Install the optional dependency when
the environment supports it:

```bash
python -m pip install ".[fftw]"
```

The runtime prints a warning and continues with NumPy if `pyFFTW` is not
available. FFmpeg is separate and is only needed when exporting video.

## Validate a change

From the repository root:

```bash
python -m pytest
python -m compileall euskera
```

For a lightweight API check:

```python
import euskera
from euskera.backgrounds import system

assert callable(euskera.evolve)
assert callable(system)
```
