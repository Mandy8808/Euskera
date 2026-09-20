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

For equations, derivations and bibliography, see
[`references/main.tex`](../references/main.tex) and
[`references/bibliografía.bib`](../references/bibliografía.bib).
The repository also contains the generated
[`SP_system.pdf`](../references/SP_system.pdf) when available. Numerical
helpers and public names are listed in the [API guide](api.md); the data flow
is described in [architecture](architecture.md).
