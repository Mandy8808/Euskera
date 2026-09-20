"""Canonical background fitting API."""
import warnings
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp, quad
from scipy.optimize import root_scalar
from euskera.observables import energy_mass as em
from .systems import systemMultifrequency
from .solvers import shoot, freq_shoot, identify, freq_shoot2, identify2

def fitting(syst, V0, indck, indXc, BCind, inddXc, limit, argf=None, info=False,
            tol=1e-14, met='RK45', Rtol=1e-07, Atol=1e-8, npt=100, klim=500):
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
    """

    # discrete points
    rmin, rmax = limit
    rspan = np.linspace(rmin, rmax, npt)
    V0 = np.array(V0, dtype='float64')

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
            ck = algebSyst(arg, info=info)

        V0[indck] = ck

        if k==klim:
            print('Stop', np.abs(Xbc-BCind), 'k', k, 'ck', V0[indck])
            break
        k += 1
    return V0

def algebSyst(arg, remNul=True, info=False):
    """
    Resuelve la ecuacion 22 de 2208.13221v1.pdf
    """
    dXc, Xbc, Xb, ck = arg

    MI = np.array(dXc, dtype='float64')
    MD = np.array(dXc@ck - (Xbc - Xb), dtype='float64')

    if remNul:  # remueve las componentes asociados al vector nulo en caso de tener
        test = MD!=0
        MI = MI[:, test][test]
        MD = MD[test]

        ck1 = np.zeros(len(test))
        temp = np.linalg.solve(MI, MD)
        ck1[test] = temp
    else:
        ck1 = np.linalg.solve(MI, MD)

    if info:
        temp = np.allclose(np.dot(MI, ck1), MD)
        print(temp)

    return ck1

__all__ = ['fitting', 'algebSyst']
