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
from euskera.evolution import EvolutionConfig
from euskera.io import OutputConfig
from euskera.observables import DiagnosticsConfig

evolution = EvolutionConfig(gridlength=10, resol=128, tmax=1)
save_config = OutputConfig(
    address="runs/my_experiment", save_number=20,
    data_save={"grid": True, "save_rho": True, "save_psi": False,
               "save_phi": False, "save_energies": True},
)
diagnostics = DiagnosticsConfig()
```

Pass `save_config` as `output_config` (and the other configurations as
`evolution_config` and `diagnostics_config`) to `euskera.evolve`. Keep the
configuration and generated data together so a run can be reconstructed.

## Select conserved quantities

The default diagnostics include particle number and energy. Additional
diagnostics can be enabled with `comp_conserv_update`:

```python
diagnostics = DiagnosticsConfig(Pi=True, Frequency=True)
```

Pass `diagnostics` as `diagnostics_config`; only enable expensive diagnostics when they are needed for the analysis,
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

The test suite also runs a short, deterministic end-to-end Gaussian evolution
with `EvolutionConfig` and `OutputConfig`, checking that output is written to
the pytest temporary directory.  It validates finite potential and density
arrays, the discrete particle-number invariant, finite RHS/profile
observables for single- and multifrequency background systems. These checks
use minimal grids and integration intervals; they are smoke tests, not
physical convergence studies. Notebook coverage is intentionally lightweight:
the suite exercises equivalent setup/API cells and never executes the
heavy, data-producing notebooks end to end.

For background diagnostics, import plotting helpers from
`euskera.visualization` and select Matplotlib's non-interactive `Agg` backend
in automated jobs.
