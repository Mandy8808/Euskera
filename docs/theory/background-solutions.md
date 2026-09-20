# Background solutions

## Scope

This page documents the radial background problem used by Euskera. The
physical equations and background families follow
[arXiv:2412.06901](https://arxiv.org/pdf/2412.06901).

The background workflow is independent from the three-dimensional dynamical
simulation workflow. A computed profile can later be used as an initial model,
but constructing the profile does not require running a simulation.

## Mathematical problem

The paper formulates radial stationary and multifrequency configurations
coupled to a gravitational potential. The implementation represents the
resulting radial systems through the functions in
[`euskera.backgrounds.systems`](../../euskera/backgrounds/systems.py).

The equations are posed with regularity conditions at the origin and
decaying or finite-energy conditions at the outer boundary. The central
amplitude, frequency parameters and polarization determine the family of
solutions. The three single-frequency polarizations are represented by the
parameters

| Polarization | `gamma` | `alpha` |
| --- | ---: | ---: |
| Linear | 0 | 0 |
| Circular | 0 | 1 |
| Radial | 1 | 0 |

The multifrequency sector uses a coupled collection of radial profiles and
the corresponding total density.

## Code correspondence

| Scientific role | Euskera module |
| --- | --- |
| Radial systems and boundary data | [`backgrounds.systems`](../../euskera/backgrounds/systems.py) |
| Frequency and node search | [`backgrounds.shooting`](../../euskera/backgrounds/shooting.py) |
| Parameter correction for a boundary-value solve | [`backgrounds.fitting`](../../euskera/backgrounds/fitting.py) |
| Profile construction and observables | [`backgrounds.profiles`](../../euskera/backgrounds/profiles.py) |
| Background mass and energy | [`observables.background_conserv_quant`](../../euskera/observables/background_conserv_quant.py) |

## References

- Chávez Nambo et al., [arXiv:2412.06901](https://arxiv.org/pdf/2412.06901),
  background solutions and their physical observables.
- The spectral stability formulation is documented separately in
  [Spectral method](spectral-method.md).

