"""MULTIFREQUENCY FITTING METHODOLOGY"""
import numpy as np
from scipy.integrate import solve_ivp

def fitting(syst, V0, indck, indXc, BCind, inddXc, limit, argf=None, info=False,
            tol=1e-14, met='RK45', Rtol=1e-07, Atol=1e-8, npt=100, klim=500,
            active=None):
    """
    syst -> System of equation with the structure: f(r, yV, arg)
            yV = [x1, x2, ..., xn, dx1c1, dx2c1, ..., dxnc1, ..., dx1cn, dx2cn, ..., dxncn]
            arg: arguments that will be passed to the system

    limit -> iteration limit: [rmin, rmax]

    V0 -> Vector with the Boundary conditions (BC).
          For example: V0 = [x1(0), c1, c2,
                             dx1c1(0), dx2c1(0), dx3c1(0),
                             dx1c2(0), dx2c2(0), dx3c2(0)]
          In this case, c1, c2 are the initial guesses or seeds for x2(0) and x3(0) respectively

    BCind -> A vector with the BC values at the right side which are used to fit the c1, ..., cn values

    indck -> A boolean vector with only True values on the c's positions on the vector V0.
            From the previous example we have:
            indck = [False, True, True,
                     False, False, False,
                     False, False, False]

    indXc -> A boolean vector with only True values on the position in yV associated to the xi variables
            with BC values at the right side.

            For example, if we have a system with the structure:
              yV = [x1, x2, x3, dx1c1, dx2c1, dx3c1, dx1c2, dx2c2, dx3c2]
            the BC:
              V0 = [x1(0), c1, c2, 0, 0, 0, 0, 0, 0]
            and:
              BCind = [x1(1), x2(1)]  notice that the x1 and x2 variables have associated a BC at the right side

            For this example:
             indXc = [True, True, False, False, False, False, False, False, False]
            the True corresponds to the respective positions in yV of the variables x1, x2

    inddXc -> A boolean vector with only True values on the position in yV associated to the derivatives of xi variables
            with BC values at the right side. From the previous example we have:
            inddXc = [False, False, False, True, True, False, True, True, False]

    argf -> arguments that will be passed to the system
    info -> Print some extra info
    active -> Optional boolean mask for active fields, in the same order as
              the selected unknowns and boundary equations. None solves the
              full system. Inactive unknowns are set to zero; their fields
              must be identically zero by the supplied initial conditions.
              This mask does not change the augmented ODE state or its masks.
    """

    # discrete points
    rmin, rmax = limit
    rspan = np.linspace(rmin, rmax, npt)
    V0 = np.array(V0, dtype='float64')
    active = _active_mask(active, int(np.sum(indck)))
    if active is not None:
        if np.sum(indXc) != active.size:
            raise ValueError("active requires one boundary equation per unknown.")
        seeds = V0[indck].copy()
        seeds[~active] = 0.0
        V0[indck] = seeds

    k = 0
    while True:
        sol = solve_ivp(syst, [rmin, rmax], V0, t_eval=rspan, args=[argf], method=met, rtol=Rtol, atol=Atol)

        Xbc = (sol.y[indXc])[:,-1]
        if np.all(np.abs(Xbc-BCind)<tol):
            print('Error', np.abs(Xbc-BCind), 'k', k, 'ck', V0[indck])
            break
        else:
            temp = (sol.y[inddXc])[:,-1]
            # notar que es necesario la traspuesta para realizar el producto de la forma adecuada: e.g. dx1*c1+dx2*c2 ...
            dXc = np.transpose(temp.reshape((np.sum(indXc), np.sum(indck))))
            arg = [dXc, Xbc, BCind, V0[indck]]
            ck = algebSyst(arg, active=active, info=info)

        V0[indck] = ck

        if k==klim:
            print('Stop', np.abs(Xbc-BCind), 'k', k, 'ck', V0[indck])
            break
        k += 1
    return V0

def _active_mask(active, size):
    if active is None:
        return None
    mask = np.asarray(active)
    if mask.dtype != np.dtype(bool) or mask.shape != (size,):
        raise ValueError("active must be a boolean mask with one entry per unknown.")
    return mask.copy()


def algebSyst(arg, active=None, info=False):
    """
    Solve the boundary correction system (equation 22 of 2208.13221v1).

    active selects corresponding rows and columns for active fields. Inactive
    unknowns are zero. None solves the full system, even when its right-hand
    side contains zeros. The mask must describe structurally inactive fields,
    not zeros caused by cancellation in the right-hand side.
    """
    dXc, Xbc, Xb, ck = arg

    MI = np.array(dXc, dtype='float64')
    ck = np.asarray(ck, dtype=float)
    if ck.ndim != 1 or MI.shape != (ck.size, ck.size):
        raise ValueError("The correction matrix must be square with one row per unknown.")
    MD = MI @ ck - (np.asarray(Xbc, dtype=float) - np.asarray(Xb, dtype=float))
    if MD.shape != ck.shape:
        raise ValueError("Boundary values must have one entry per unknown.")
    active = _active_mask(active, ck.size)

    if active is not None:
        ck1 = np.zeros(ck.size)
        if np.any(active):
            ck1[active] = np.linalg.solve(MI[np.ix_(active, active)], MD[active])
    else:
        ck1 = np.linalg.solve(MI, MD)

    if info:
        temp = np.allclose(np.dot(MI, ck1), MD)
        print(temp)

    return ck1

__all__ = ['fitting', 'algebSyst']
