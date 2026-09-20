"""Canonical background systems API."""
import warnings
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp, quad
from scipy.optimize import root_scalar
from euskera.observables import energy_mass as em

def system(r, yV, arg):
    r"""
    System of equations (52) for one scalar field. Note that we used the Ansatz:
        \sigma = r^\gamma \sigma, with \gamma = 0 (multi-frequency, Lineal, Circular), \gamma = 1 (radial)
    in order to incorporated the radial case.

    Variables:
        yV = [p0, p1, u0, u1]  # System components
        arg = [sumpComp, LambT, gamma]  # System parameters

    Returns:
        [f0, f1, f2, f3]
    """

    # Ensure yV has exactly four elements
    if len(yV) != 4:
        raise ValueError("yV must contain exactly four elements: [p0, p1, u0, u1]")

    p0, p1, u0, u1 = yV
    sumpComp, LambT, gamma = arg

    # Validate parameters
    if sumpComp < 0:
        raise ValueError("sumpComp cannot be negative.")

    # If p0 is too large, return zeros to avoid instability
    if np.abs(p0) > 80:
        return np.array([0, 0, 0, 0])

    # Evaluate the system based on r
    if r > 0:
        f0 = p1
        f1 = LambT*sumpComp*p0*r**(2*gamma) - 2*(1 + gamma)*p1/r - u0*p0
        f2 = u1
        f3 = -r**(2*gamma)*sumpComp - 2*u1/r
    elif r == 0:
        f0 = p1
        f1 = (LambT*sumpComp*p0*r**(2*gamma) - u0*p0)/(2*gamma + 3)
        f2 = u1
        f3 = -r**(2*gamma)*sumpComp/3
    else:
        raise ValueError("r cannot be negative in this context.")

    return [f0, f1, f2, f3]

def systemMultifrequency(r, yV, arg):
    r"""
    Full System of equations (52), with the Ansatz
    \sigma = r^\gamma \sigma, with \gamma = 0 (multi-frequency, Lineal, Circular), \gamma = 1 (radial)
    in order to incorporated the radial case.

    Parameters:
        r   : float - Radial coordinate
        yV  : array - System components of size [4 * numFields]
        arg : tuple - (numFields, LambT, gamma)

    Returns:
        Array of derivatives for each component
    """

    if len(arg) == 1:
        gamma = 0
        numFields = int(len(yV)/4)
        LambT, = arg
    elif len(arg) == 2:
        numFields = int(len(yV)/4)
        LambT, gamma = arg
    else:
        numFields, LambT, gamma = arg

    # Checking that for multifrequency case, gamma must be zero
    if numFields != 1 and gamma != 0:
        raise ValueError("gamma must be exactly zero for the multifrequency case.")

    # Ensuring yV has the correct number of elements
    if len(yV) != 4 * numFields:
        raise ValueError("yV must have exactly %d elements." % (4 * numFields))

    # Compute the sum of squared components
    sumpComp = np.sum(yV[:2*numFields:2]**2)

    # Preallocate result array
    valores = np.zeros(4 * numFields)

    # Handling extremely large values to prevent numerical instability
    if sumpComp > 1e6:
        warnings.warn("The density is extremely large and will cause numerical instabilities.")
        # raise ValueError("The density is extremely large and will cause numerical instabilities.")
        return valores

    # Compute derivatives for each component
    arg_Componente = [sumpComp, LambT, gamma]
    for componente in range(0, 2*numFields, 2):
        p0, p1 = yV[componente], yV[componente+1]
        u0, u1 = yV[2*numFields + componente], yV[2*numFields + componente+1]
        yV_Componente = [p0, p1, u0, u1]
        f0, f1, f2, f3 = system(r, yV_Componente, arg_Componente)
        valores[componente], valores[componente+1] = f0, f1
        valores[2*numFields + componente], valores[2*numFields + componente+1] = f2, f3

    return valores

def systemMultFreqTot(r, yVT, arg):
    """
    syst -> System of equation with the structure: f(r, yV, arg)
            yV = [x1, x2, ..., xn, dx1c1, dx2c1, ..., dxnc1, ..., dx1cn, dx2cn, ..., dxncn]
            arg: arguments that will be passed to the system

    argf -> arguments that will be passed to the system
    info -> Print some extra info

    yVT -> [p0x, p1x, p0y, p1y, p0z, p1z, u0x, u1x, u0y, u1y, u0z, u1z,\
           Duxp0x, Duxp1x, Duxp0y, Duxp1y, Duxp0z, Duxp1z, Duxu0x, Duxu1x, Duxu0y, Duxu1y, Duxu0z, Duxu1z,\
           Duyp0x, Duyp1x, Duyp0y, Duyp1y, Duyp0z, Duyp1z, Duyu0x, Duyu1x, Duyu0y, Duyu1y, Duyu0z, Duyu1z,\
           Duzp0x, Duzp1x, Duzp0y, Duzp1y, Duzp0z, Duzp1z, Duzu0x, Duzu1x, Duzu0y, Duzu1y, Duzu0z, Duzu1z]
    """

    yV = yVT[:12]
    f0to11 = systemMultifrequency(r, yV, arg)   # [f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11]

    # Derivadas con respecto a u_i
    Mat = MatrizDXDu(r, yV, arg, info=False)
    yVDux = yVT[12:24]
    yVDuy = yVT[24:36]
    yVDuz = yVT[36:]

    f12tof23 = Mat@yVDux  # Dux [f12, f13, f14, f15, f16, f17, f18, f19, f20, f21, f22, f23]
    f24to35 = Mat@yVDuy  # Duy [f24, f25, f26, f27, f28, f29, f30, f31, f32, f33, f34, f35]
    f36to47 = Mat@yVDuz  # Duz  [f36, f37, f38, f39, f40, f41, f42, f43, f44, f45, f46, f47]

    yVOut = np.concatenate((f0to11, f12tof23, f24to35, f36to47))

    return yVOut

def MatrizDXDu(r, yV, arg, info=False):
    r"""
    Matriz creada a partir de computar d/dz(\partial X/\partial c_i)
    ver 2208.13221v1.pdf
    """
    #numFields, LambT, gamma = arg
    LambT, = arg
    p0x, p1x, p0y, p1y, p0z, p1z, u0x, u1x, u0y, u1y, u0z, u1z = yV

    sumpi = p0x**2 + p0y**2 + p0z**2

    # llenando matriz
    if r==0:
        Mat = np.zeros((12, 12))
    else:
        Mat = np.diag([0, -2/r]*6)  # lleno la diagonal

    diag1S = np.diag([1, 0]*6, k=1)[:-1,:-1]
    Mat += diag1S

    temp1 = np.array([-p0x, -p0y, -p0z])
    Mat[7::2, :5:2] = 2*temp1
    Mat[1, 6] = temp1[0]
    Mat[3, 8] = temp1[1]
    Mat[5, 10] = temp1[2]

    t1, t2, t3 = 2*LambT*p0y*p0x, 2*LambT*p0z*p0x, 2*LambT*p0z*p0y
    temp2 = [
        [LambT*(sumpi+2*p0x**2)-u0x, t1, t2],
        [t1, LambT*(sumpi+2*p0y**2)-u0y, t3],
        [t2, t3, LambT*(sumpi+2*p0z**2)-u0z]]
    Mat[1:6:2, :5:2] = temp2

    if r==0:
        #np.fill_diagonal(Mat, 0)
        Mat = Mat/3

    if info:
        print(Mat)

    return Mat

__all__ = ['system', 'systemMultifrequency', 'systemMultFreqTot', 'MatrizDXDu']
