# Architecture and simulation flow

## Package boundaries

Euskera is one scientific package with two primary workflows:

| Namespace | Responsibility |
| --- | --- |
| `euskera.evolution` | Three-dimensional field evolution and potential solving |
| `euskera.backgrounds` | Radial background solutions, shooting and profile construction |

Shared concerns have explicit namespaces:

| Namespace | Responsibility |
| --- | --- |
| `euskera.core` | Grids and shared field primitives |
| `euskera.models` | Initial field models |
| `euskera.numerics` | General numerical methods |
| `euskera.spectral` | Spectral discretisation, operators, and eigenvalue methods |
| `euskera.observables` | Energy, mass and conserved quantities |
| `euskera.io` | Simulation input/output |
| `euskera.visualization` | Plots and videos |

There are no compatibility packages outside these canonical namespaces;
notebooks and downstream users should import the workflow boundaries above.

## Evolution pipeline

`euskera.evolve` (also available as `euskera.evolution.evolve`) coordinates
the following stages:

1. **Model selection.** The `Models` registry resolves one or more model names
   and validates their parameter dictionaries.
2. **Grid construction.** Real-space and Fourier-space grids are built from
   `gridlength` and `resol`.
3. **Initial field.** The selected model creates the wavefunction `psi` and
   initial density `rho_i`.
4. **Initial potential.** `Upotential` solves the Poisson contribution in
   Fourier space.
5. **Time evolution.** `PKP` advances the fields with the split-step
   Fourier method.
6. **Diagnostics and output.** Conserved quantities are computed according to
   `comp_conserv`; snapshots and parameters are written using the saving
   configuration.

The flow is intentionally composable: model construction, grid generation and
diagnostic functions can be called independently in a notebook or test.

## Model registry

The built-in model names are:

| Name | Use |
| --- | --- |
| `soliton` | Radial soliton profiles |
| `gaussian_function` | One or more Gaussian field contributions |
| `ell_boson` | Ellipsoidal boson-star configurations |
| `proca` | Proca-star configurations and polarisation data |

Models are applied in the order in which they are supplied to `Models`. This
allows a base configuration and additional perturbations to be composed.

## Output flow

By default, simulations save NPZ data under the configured `address` (the
default is `Data`) and write a `parameters` file describing the run. The
`data_save` mapping controls whether full fields, planes, lines and energies
are saved. HDF5 is available through the saving API when `h5py` is installed.

Avoid writing generated output into the source directories tracked by Git.
Use a dedicated local output directory and retain the parameter file with the
resulting data.

Background plotting helpers are canonical visualization exports rather than
modules under `euskera.backgrounds`; shared Matplotlib configuration is also
centralized in `euskera.visualization.configuration`.
