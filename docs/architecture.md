# Architecture and simulation flow

## Package boundaries

Euskera is distributed as two importable packages:

| Package | Responsibility |
| --- | --- |
| `euskera` | Grid construction, initial models, time evolution, conserved quantities, plotting and data storage |
| `background` | Radial background solutions, spectral methods, energy/mass calculations and background plotting |

The top-level `euskera` module re-exports the most commonly used symbols.
Imports from submodules are useful when a workflow needs a more focused
dependency, for example `from euskera.main.grids import RealGrid`.

## Evolution pipeline

`euskera.evolve` coordinates the following stages:

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
