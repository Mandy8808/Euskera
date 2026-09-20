# Background shooting workflow

This workflow constructs radial background solutions from the equations and
boundary conditions described in [Background solutions](../theory/background-solutions.md).

```mermaid
flowchart TD
    Parameters[Physical parameters]
    System[Radial system]
    Guess[Frequency or initial guess]
    Integrate[Radial integration]
    Events[Node and event detection]
    Update[Bracket or midpoint update]
    Profile[Background profile]

    Parameters --> System
    System --> Guess
    Guess --> Integrate
    Integrate --> Events
    Events --> Update
    Update --> Integrate
    Events --> Profile
```

The shooting stage and the fitting stage are separate. Shooting selects the
frequency or interval that produces the desired event and node structure.
Fitting, when needed, corrects unknown initial parameters using the
boundary-value method described in
[Boundary-value numerical method](../methods/boundary-value-method.md).

Relevant modules:

```python
from euskera.backgrounds import profiles, shooting, systems
```

After a profile is constructed, use the background observables and plotting
helpers to inspect it before passing it to a dynamical simulation.

