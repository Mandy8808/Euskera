# Spectral method and stability

## Scope

The spectral formulation used for radial operators and perturbation spectra
follows [arXiv:2512.04376](https://arxiv.org/pdf/2512.04376). It is not the
same operation as constructing a background profile with shooting or fitting.

The spectral workflow takes a background as input, discretizes radial
operators with Chebyshev points, builds the perturbation blocks and solves an
eigenvalue problem.

## Computational flow

```mermaid
flowchart LR
    Profile[Background profile]
    Grid[Chebyshev grid]
    Operators[Radial operators]
    Blocks[Polarization blocks]
    Spectrum[Eigenvalues and modes]

    Profile --> Grid
    Grid --> Operators
    Operators --> Blocks
    Blocks --> Spectrum
```

## Code correspondence

| Scientific role | Euskera module |
| --- | --- |
| Chebyshev points and differentiation matrices | [`spectral.chebyshev`](../../euskera/spectral/chebyshev.py) |
| Background radial operators | [`spectral.operators`](../../euskera/spectral/operators.py) |
| Polarization and multifrequency blocks | [`spectral.blocks`](../../euskera/spectral/blocks.py) |
| Eigenvalue assembly and solution | [`spectral.eigensolver`](../../euskera/spectral/eigensolver.py) |
| Eigenvalue organization | [`spectral.analysis`](../../euskera/spectral/analysis.py) |
| Spectral plots | [`visualization.spectral_plot`](../../euskera/visualization/spectral_plot.py) |

The public numerical entry points are:

```python
from euskera.spectral import LamJval, backgroundOper, cheb, spectrum
```

## Interpretation

The eigenvalues describe linear perturbation modes around a previously
constructed background. The spectral calculation therefore belongs after the
background workflow and before any interpretation of stability. It does not
replace the radial solver and it does not fit the background data.

