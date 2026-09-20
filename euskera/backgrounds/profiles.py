"""BACKGROUND PROFILES API"""
import numpy as np
from scipy.integrate import solve_ivp
from euskera.observables import energy_mass as em
from .systems import systemMultifrequency

def profilesFromSolut(datos, mult=False, rmin=0, fac=4*np.pi, Nptos=2000, info=False,
                      met='DOP853', Rtol=1e-09, Atol=1e-10):  # 'RK45'
    """
    Generates discrete profiles from a numerical solution.

    Parameters:
    -----------
    datos : tuple
        Solution data, with different structures depending on the `mult` value.
    mult : bool, optional
        Indicates whether it involves multiple fields (default `False`).
    rmin : float, optional
        Minimum integration radius (default `0`).
    Nptos : int, optional
        Number of points in the discretization (default `2000`).
    info : bool, optional
        If `True`, prints the mass and energy values.

    Returns:
    --------
    dict with the following keys:
        - "energy": Computed energy.
        - "mass": Computed mass.
        - "profiles_f": Profiles of `f`.
        - "profiles_df": Derivatives of `f`.
        - "profiles_u": Profiles of `u`.
        - "profiles_du": Derivatives of `u`.
        - "LambT": LambdaT parameter.
        - "gamma": Gamma parameter.
    """

    # Extracting data
    if mult:
        LambT, nodes, f0 = datos[-1], datos[-2], datos[-3]
        numFields = len(nodes)
        u0, rmax = datos[:numFields], datos[numFields]
        gamma = 0
    else:
        f0, rmax, gamma, LambT, nodes, _, met, Rtol, Atol, u0 = datos
        numFields = 1

    # Initial conditions setup
    V0 = np.zeros(4 * numFields)
    V0[:2 * numFields:2] = f0
    V0[2 * numFields::2] = u0

    rspan = np.linspace(rmin, rmax, Nptos)
    arg = [numFields, LambT, gamma]

    sol = solve_ivp(systemMultifrequency, [rmin, rmax], V0, t_eval=rspan,
                     args=[arg], method=met, rtol=Rtol, atol=Atol)

    # Extracting profiles
    perfSig = [sol.y[i] for i in range(0, 2 * numFields, 2)]
    perfdSig = [sol.y[i] for i in range(1, 2 * numFields, 2)]
    perfU = [sol.y[i] for i in range(2 * numFields, 4 * numFields, 2)]
    perfdU = [sol.y[i] for i in range(2 * numFields + 1, 4 * numFields, 2)]

    # Computing energy and mass
    rD = sol.t
    sigtot = sum([perfSig[i]**2 for i in range(numFields)])
    energy, mass = [], []
    for V0, sig in zip(V0[2 * numFields::2], perfSig):
        En = em.energEng(rD, sigtot, V0, gamma=gamma, kind='quadratic', fill_value="extrapolate")
        Mas = em.massVal(rD, sigtot, gamma=gamma, fac=fac, kind='quadratic', fill_value="extrapolate")
        energy.append(En)
        mass.append(Mas)

    if info:
        print(f"Mass values: {mass}")
        print(f"Eigenvalues of energy: {energy}")

    return energy, mass, rD, perfSig, perfdSig, perfU, perfdU, LambT, gamma

__all__ = ['profilesFromSolut']
