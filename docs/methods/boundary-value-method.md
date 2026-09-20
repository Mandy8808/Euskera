# Boundary-value numerical method

## Reference

The iterative boundary-value methodology implemented by the background
helpers is described in [arXiv:2208.13221](https://arxiv.org/pdf/2208.13221).
That paper discusses shooting formulations, Newton-type corrections and the
linear algebra used to enforce conditions at the far boundary.

This reference is a numerical-method reference. It does not define the
physical background equations; those are documented in
[Background solutions](../theory/background-solutions.md).

## Two distinct operations

Euskera keeps two operations separate:

### Shooting

[`backgrounds.shooting`](../../euskera/backgrounds/shooting.py) searches
frequency intervals and identifies the requested node structure from event
crossings. Its correction logic uses bracketing and midpoint updates.

### Fitting

[`backgrounds.fitting`](../../euskera/backgrounds/fitting.py) integrates an
augmented system containing the original variables and sensitivity
derivatives. It evaluates the residual at the right boundary and calls
`algebSyst` to update the unknown initial parameters.

This fitting function is a boundary-condition parameter-correction step. It
must not be confused with fitting a closed-form curve to observational data.

## Fitting iteration

The implementation follows this structure:

```text
Initial values and parameter guesses
                ↓
Integrate the augmented ODE system
                ↓
Evaluate the right-boundary residual
                ↓
Build the correction system
                ↓
Update the unknown initial parameters
                ↓
Repeat until tolerance or iteration limit
```

The main controls are:

- `tol`: boundary residual tolerance;
- `Rtol`, `Atol`: ODE integration tolerances;
- `npt`: number of evaluation points;
- `klim`: maximum correction iterations;
- `met`: `solve_ivp` integration method.

The method returns the corrected initial vector `V0`, which can then be used
by the profile workflow.

