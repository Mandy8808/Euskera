"""Multifrequency background shooting workflows."""
import warnings
import numpy as np
from scipy.integrate import solve_ivp
from euskera.visualization import background_plots as cplot
from .systems import systemMultifrequency
from .solvers import shoot, freq_shoot, identify, freq_shoot2, identify2

def MultFreq_solveG2_vO(Ini0, Uintrs, rmax, rmin=0, LambT=1, nodos=[0, 0, 0],
                     met='RK45', Rtol=1e-09, Atol=1e-10, lim=1e-6, info=False,
                     klim=500, outval=13, delta=0.4): #'DOP853''LSODA'
    """
    In:
    Uintx -> [Umin, Umax]
    Uinty -> [Umin, Umax]
    Uintz -> [Umin, Umax]
    Ini0 -> [p0x, p1x, p0y, p1y, p0z, p1z, u1x, u1y, u1z]
    rmax, rmin ->
    LambT -> +1 Repulsive case, -1 Atractive case
    nodos -> [nodos_p0x, nodos_p0y, nodos_p0z]

    Orden de las variables
    [phix, phix', phiy, phiy', phiz, phiz', ux, ux', uy, uy', uz, uz'] -> [p0x, p1x, p0y, p1y, p0z, p1z, u0x, u1x, u0y, u1y, u0z, u1z]
    """

    nodos = np.array(nodos)
    p0x, p1x, p0y, p1y, p0z, p1z, u1x, u1y, u1z = Ini0

    p0Data = [p0x, p0y, p0z]
    Uintx, Uinty, Uintz = Uintrs
    Uminx, Umaxx = Uintx
    Uminy, Umaxy = Uinty
    Uminz, Umaxz = Uintz

    print('Finding a profile with nx, ny, nz', nodos, 'nodes')

    # Events
    def Sigx(r, U, arg): return U[0]
    def dSigx(r, U, arg): return U[1]
    def Sigy(r, U, arg): return U[2]
    def dSigy(r, U, arg): return U[3]
    def Sigz(r, U, arg): return U[4]
    def dSigz(r, U, arg): return U[5]
    Sigx.direction = 0; dSigx.direction = 0
    Sigy.direction = 0; dSigy.direction = 0
    Sigz.direction = 0; dSigz.direction = 0

    # ordenando de mayor a menor para la iteracion
    k = 0

    arg = [LambT]
    Uintrs = np.array([[Uminx, Umaxx], [Uminy, Umaxy], [Uminz, Umaxz]])
    sigModifold = None
    UintrOrig = np.copy(Uintrs)
    out = 0
    while True:
        u0 = np.array([shoot(*i) for i in Uintrs])
        V0 = [p0x, p1x, p0y, p1y, p0z, p1z, u0[0], u1x, u0[1], u1y, u0[2], u1z]

        sol = solve_ivp(systemMultifrequency, [rmin, rmax], V0, events=(Sigx, dSigx, Sigy, dSigy, Sigz, dSigz),
                         args=(arg,), method=met,  rtol=Rtol, atol=Atol)

        eventos = np.array([[sol.t_events[0], sol.t_events[1]],
                   [sol.t_events[2], sol.t_events[3]],
                   [sol.t_events[4], sol.t_events[5]]], dtype=object)
        sigModif = identify2(sol.t_events, nodos, p0Data, info=info)

        if info:
            print(sigModif)

        iInterv, rTemp = freq_shoot2(eventos[sigModif], nodos[sigModif], u0[sigModif], Uintrs[sigModif], rmax)

        if abs((iInterv[1]-iInterv[0])/2) <= lim:
            if info:
                print(out)
                print('Maxima precisión alcanzada: U0x = ', V0[6], ' U0y = ', V0[8], ' U0z = ', V0[10], 'radio', rTemp)

            if out==outval:
                print('Maxima precisión alcanzada: U0x = ', V0[6], ' U0y = ', V0[8], ' U0z = ', V0[10], 'radio', rTemp)
                u0 = [V0[6], V0[8], V0[10]]
                return u0, rTemp, sol.t_events[::2]
            else:
                #Uintrs = np.copy(UintrOrig) # reinicio los que  ya no son iguales
                Uintrs[sigModif] = [iInterv[0]-delta, iInterv[1]+delta]
                out += 1
        else:
            Uintrs[sigModif] = iInterv

        if info:
            print(Uintrs)

        if np.all(np.array([shoot(*i) for i in Uintrs])==u0):
            print('Found: U0x = ', V0[6], ' U0y = ', V0[8], ' U0z = ', V0[10], 'radio', rTemp)
            u0 = [V0[6], V0[8], V0[10]]
            return u0, rTemp, sol.t_events[::2]

        if k==klim:
            print('loop limit reached')
            break

        k += 1

def MultFreq_solveG2(Ini0, numFields, Uintrs, rmax, rmin=0, LambT=0, gamma=0, nodos=None,
                     met='RK45', Rtol=1e-09, Atol=1e-10, lim=1e-14, info=False,
                     klim=500, outval=13, delta=0.4, part=None, inv=False):
    """
    Multi-frequency solver using solve_ivp.
    """

    # Ensuring gamma is zero for multifrequency case
    if numFields != 1 and gamma != 0:
        raise ValueError("gamma must be exactly zero for the multifrequency case.")

    # Validating initial conditions
    if len(Ini0) != 3 * numFields:
        raise ValueError(f"Ini0 must have exactly {3 * numFields} elements.")

    if len(Uintrs) != numFields:
        raise ValueError(f"Uintrs must have exactly {numFields} intervals.")

    if sum([len(Uintrs[i]) == 2 for i in range(numFields)]) != numFields:
        raise ValueError("Each Uintrs component must have exactly 2 extremal values.")


    Uintrs = np.array(Uintrs, dtype=object)  # Ensure it's modifiable
    arg = [numFields, LambT, gamma]

    # Setting node values
    if nodos is None:
        nodos = np.zeros(numFields, dtype=np.int8)
    else:
        if len(nodos) != numFields:
            raise ValueError(f"The node number must have exactly {numFields} elements.")
        nodos = np.array(nodos, dtype=np.int8)

    print('Finding a profile with nodes:', nodos)

    # Define event functions based on numFields
    def Sigx(r, U, arg): return U[0]
    def dSigx(r, U, arg): return U[1]

    Sigx.terminal = False
    dSigx.terminal = False

    FindEvents = [Sigx, dSigx]

    if numFields == 2:
        def Sigy(r, U, arg): return U[2]
        def dSigy(r, U, arg): return U[3]

        Sigy.terminal = False
        dSigy.terminal = False
        FindEvents.extend([Sigy, dSigy])

    if numFields == 3:
        def Sigz(r, U, arg): return U[4]
        def dSigz(r, U, arg): return U[5]

        Sigz.terminal = False
        dSigz.terminal = False
        FindEvents.extend([Sigz, dSigz])

    # Initial conditions setup
    V0 = np.zeros(4 * numFields)
    V0[:2 * numFields] = Ini0[:2 * numFields]
    V0[2 * numFields + 1::2] = Ini0[2 * numFields:]

    # Track the central amplitude values
    p0Data = [Ini0[i] for i in range(0, 2 * numFields, 2)]

    k = 0
    out = 0
    while k < klim:
        u0 = np.array([shoot(*i) for i in Uintrs])
        V0[2 * numFields::2] = u0  # Updating boundary conditions

        # Solving ODE
        sol = solve_ivp(systemMultifrequency, [rmin, rmax], V0, events=FindEvents,
                        args=(arg,), method=met, rtol=Rtol, atol=Atol)

        if info:
            cplot.plotUsingSol(sol)

        # Extracting event times
        eventos = np.array([sol.t_events[i:i+2] for i in range(0, 2*numFields, 2)], dtype=object)

        if part:
            even, shot = part
            iInterv, rTemp = freq_shoot(eventos[even], nodos[even], u0[shot], Uintrs[shot], rmax)
            sigModif = shot
        else:
            sigModif_bool = identify(sol.t_events, nodos, p0Data)
            sigModif = np.argmax(sigModif_bool)  # Convert boolean list to index

            iInterv, rTemp = freq_shoot(eventos[sigModif], nodos[sigModif], u0[sigModif], Uintrs[sigModif], rmax, inv=inv)

        if info:
            print("Modifying:", sigModif)


        if abs((iInterv[1] - iInterv[0]) / 2) <= lim:
            if info:
                print(out, '->>', iInterv[1], iInterv[0])
                print('Precision reached: U0 =', u0, 'radius', rTemp)

            if out == outval:
                print('Final precision: U0 =', u0, 'radius', rTemp)
                return u0, rTemp, sol.t_events[::2]
            else:
                Uintrs[sigModif] = [iInterv[0] - delta, iInterv[1] + delta]
                out += 1
        else:
            Uintrs[sigModif] = iInterv

        if info:
            print("Updated Uintrs:", Uintrs)

        if np.all(np.array([shoot(*i) for i in Uintrs]) == u0):
            print('Solution found: U0 =', u0, 'radius', rTemp)
            return u0, rTemp, sol.t_events[::2]

        k += 1

    print('Loop limit reached')
    return None

__all__ = ['MultFreq_solveG2_vO', 'MultFreq_solveG2']
