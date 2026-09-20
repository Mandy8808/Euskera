# Scientific background

Euskera evolves self-gravitating fields on a three-dimensional Cartesian
grid. The Schrödinger–Poisson (SP) system describes a non-relativistic
self-gravitating scalar field; the Gross–Pitaevskii–Poisson (GPP) extension
adds nonlinear self-interaction. The evolution uses a split-step Fourier
method and solves the potential spectrally.

Initial configurations include solitons, Gaussian perturbations, ellipsoidal
boson stars and Proca stars (including multi-frequency cases). The separate
`euskera.backgrounds` workflow provides radial shooting, spectral profiles
and background diagnostics.

The background and simulation workflows are independent:

1. [Background solutions](theory/background-solutions.md) constructs radial
   profiles from the equations and boundary conditions.
2. [Dynamical simulation](workflows/simulation.md) evolves an initial field
   in three dimensions.

Within the background workflow, [shooting](workflows/background-shooting.md)
selects frequencies and node structures, while
[fitting](workflows/background-fitting.md) performs a separate
boundary-condition parameter correction.

The main scientific references are:

- [arXiv:2412.06901](https://arxiv.org/pdf/2412.06901) for background
  solutions and observables;
- [arXiv:2512.04376](https://arxiv.org/pdf/2512.04376) for the spectral
  perturbation and stability method;
- [arXiv:2208.13221](https://arxiv.org/pdf/2208.13221) for the numerical
  boundary-value methodology implemented by the fitting helper.

For additional equations, derivations and bibliography, see
[`references/main.tex`](../references/main.tex) and
[`references/bibliografía.bib`](../references/bibliografía.bib).
The repository also contains the generated
[`SP_system.pdf`](../references/SP_system.pdf) when available. Numerical
helpers and public names are listed in the [API guide](api.md); the data flow
is described in [architecture](architecture.md).

Background solvers do not own plotting configuration. For diagnostic plots,
use `euskera.visualization.background_plots`; this keeps solver modules
independent of plotting implementation details.
