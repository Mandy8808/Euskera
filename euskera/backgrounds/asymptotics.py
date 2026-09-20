"""Canonical background asymptotics API."""
import warnings
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp, quad
from scipy.optimize import root_scalar
from euskera.observables import energy_mass as em
from . import plot_conf as cplot
from .profiles import profilesFromSolut

def mainExt(datos, mult=False, rmin=0, Nptos=2000, fac=1, info=False,
            met='DOP853', metExp=1, Rtol=1e-09, Atol=1e-10, deltaExt=100, indF=-160,
            xscale='linear', yscale='linear', sigfilt=1e-28):
    """
    Function to numerically solve and extend profiles.

    IMPORTANT: fac=4*np.pi give the "True" mass, but for the extention we need to used fac=1
    """

    # Solving numerically using datos
    En, Mas, rD, perfSig, perfdSig, perfU, perfdU, _, gamma = profilesFromSolut(
        datos, mult=mult, rmin=rmin, Nptos=Nptos, met=met, Rtol=Rtol, Atol=Atol, fac=fac, info=info)

    numFields = len(perfSig)

    # filtering data
    sigtot = sum([perfSig[i]**2 for i in range(numFields)])
    ind = sigtot/sigtot[0] > sigfilt
    rD = rD[ind]
    perfSig = [Sig[ind] for Sig in perfSig]
    perfdSig = [dSig[ind] for dSig in perfdSig]
    perfU = [fU[ind] for fU in perfU]
    perfdU = [dfU[ind] for dfU in perfdU]

    # Making a datosOrg matrix array
    datosOrg = np.empty((numFields, 6), dtype=object)  # Using object dtype to hold lists/arrays

    temp = [En, Mas, perfSig, perfdSig, perfU, perfdU]
    for i in range(6):
        datosOrg[:, i] = temp[i]

    # Extending profiles
    dataEx = []
    for i in range(numFields):
        Ei, Mi, sD, dsD, uD, duD = datosOrg[i]

        if np.allclose(sD, 0):  # Check if sD is all zeros
            continue  # Avoid extending a null profile

        Ext = (sD[0] * rD[-1]) + deltaExt
        Np = int(Ext / 2)

        rDnew, sDnew, dsDnew, uDnew, duDnew = extend(
            gamma=gamma, rD=rD, sD=sD, dsD=dsD, uD=uD, duD=duD,
            Ext=Ext, En=Ei, Mas=Mi, met=metExp, Np=Np, inf=info, indF=indF
        )

        if info:
            # Plotting profiles for debugging/visualization
            perf1 = [rDnew, sDnew, dsDnew, uDnew, duDnew]
            perf2 = [rD, sD, dsD, uD, duD]
            cplot.plotUsingDiscSol(perf1, perf2=perf2, xscale=xscale, yscale=yscale)

        dataEx.append([Ei, rDnew, sDnew, dsDnew, uDnew, duDnew])

    return dataEx

def extend(gamma, rD, sD, dsD, uD, duD, Ext, En=None, Mas=None,
           met=1, Np=1000, inf=False, fac=1, indF=-160):
    """
    Extend radial data by an additional segment using one of two methods.

    Parameters:
    - rD: Radial data (initial).
    - sD: Profile data for the initial segment.
    - dsD: Derivative of profile data for the initial segment.
    - uD: Another data set to extend.
    - duD: Derivative of uD.
    - Ext: Length to extend the radial data.
    - En: Energy, required only for `met == 2` (optional).
    - Mas: Mass, required only for `met == 2` (optional).
    - met: Method selector (1 for Met1, 2 for Met2).
    - Np: Number of points for the extension.
    - inf: Information flag to control verbosity (for Met2).

    Returns:
    - rDnew: Extended radial data.
    - sDnew: Extended profile data.
    - dsDnew: Extended derivative of profile data.
    - uDnew: Extended data.
    - duDnew: Extended derivative of uD.

    IMPORTANT: fac=4*np.pi give the "True" mass, but for the extention we need to used fac=1
    """

    # Ensure `En` and `Mas` are provided when using method 2
    if En is None or Mas is None:
        warnings.warn("En and Mas are not provided. We computed them from the field profile.\n IMPORTANT: for multifrequency, this procedure is not correct.")

        if gamma not in [0, 1]:
            raise ValueError("Invalid value for `gamma`. It must be 0 or 1.")

        sigtot = sD**2
        V0 = uD[0]

        # Compute `En` and `Mas`
        En = em.energEng(rD, sigtot, V0, gamma=gamma, kind='quadratic', fill_value="extrapolate")
        Mas = em.massVal(rD, sigtot, gamma=gamma, fac=fac, kind='quadratic', fill_value="extrapolate")

    # Select method based on `met`
    if met == 2:
        rDnew, sDnew, dsDnew, uDnew, duDnew = extend_Met2(En, Mas, rD[:indF], sD[:indF], dsD[:indF], uD[:indF], duD[:indF], Ext,
                                                           Np=Np, info=inf)
    else:
        rDnew, sDnew, dsDnew, uDnew, duDnew = extend_Met1(En, Mas, rD[:indF], sD[:indF], dsD[:indF], uD[:indF], duD[:indF], Ext,
                                                           Np=Np)

    return rDnew, sDnew, dsDnew, uDnew, duDnew

def sParam(r, S, B=0):
    """
    Computes parameters C, k, and s based on the last two elements of r and S.

    Parameters:
        r (array-like): A sequence of radial values.
        S (array-like): A sequence of corresponding function values.
        B (float): Mass of the configuration. When B=0 the result correspond to met=1, else met=2

    Returns:
        tuple: (C, k, s), where
            - C is a scaling constant,
            - k is a decay/growth rate
    """
    # Extracting the last two elements
    r1, r2 = r[-2], r[-1]
    yr1, yr2 = S[-2], S[-1]

    # Avoid division by zero
    if yr2 == 0 or r2 == 0:
        raise ValueError("Division by zero encountered in logarithm computation.")

    # Compute k and C
    ratio = yr1 * r1 / (yr2 * r2)

    if ratio <= 0:
        #raise ValueError("Logarithm of a non-positive number encountered.")
        print(ValueError("Logarithm of a non-positive number encountered."))

    k = np.real(np.log(np.abs(ratio)))

    # Notice that if B=0 the result correspond to the met=1 else met=2
    s = np.exp(-k * r1)/r1**(1 - B/(2 * k))
    C = yr1/s  # C = yr1/s

    return C, k

def Asy_uProf(r, A, B):
    """
    Computes the asymptotic potential profile: U = A + B/r and its derivative dy.

    Parameters:
        r (float or np.ndarray): Radial coordinate (must be nonzero).
        A (float): Constant term (eingenvalue of the energy).
        B (float): Mass of the configuration.

    Returns:
        tuple: (y, dy), where
            - y  = A + B/r
            - dy = derivative of y with respect to r
    """
    # Avoid division by zero
    if np.any(r == 0):
        raise ValueError("r must be nonzero to avoid division by zero.")

    # Compute y and dy
    y = A + B/r
    dy = -B/r**2

    return y, dy

def Asy_sProf(r, C, k):
    """
    Computes the asymptotic field profile: sigma = C*exp(-k*r)/r and its derivative dy.

    Parameters:
        r (float or np.ndarray): Radial coordinate (must be nonzero).
        C (float): Scaling constant.
        k (float): Decay/growth rate.

    Returns:
        tuple: (y, dy), where
            - y  = C * exp(-k*r) / r
            - dy = derivative of y with respect to r
    """
    # Avoid division by zero
    if np.any(r == 0):
        raise ValueError("r must be nonzero to avoid division by zero.")

    # Compute y and dy
    exp_term = np.exp(-k * r)
    y = C * exp_term / r
    dy = -(C * exp_term * (1 + k * r)) / r**2

    return y, dy

def Asy_sProf_v2(r, C, En, M):
    """
    Computes a modified decaying exponential function with an additional power-law factor.

    Parameters:
    - r: float or np.array, radial coordinate (must be > 0 to avoid division errors)
    - C: float, amplitude constant
    - En: float, energy parameter (must be <= 0 to avoid sqrt of negative number)
    - M: float, mass parameter

    Returns:
    - y: computed function value
    - dy: first derivative
    """
    if En > 0:
        raise ValueError("En must be non-positive (≤ 0) to avoid complex values.")

    if np.any(r <= 0):
        raise ValueError("r must be strictly positive to prevent division errors.")

    # Compute y and dy
    k = np.sqrt(-En)
    exp_factor = np.exp(-k * r)
    power_factor = np.power(r, 1 - M/(2 * k))
    power_factor_derivada = np.power(r, -2 + M/(2 * k))

    y = C * exp_factor / power_factor
    dy = C * exp_factor * power_factor_derivada * (M - 2*k*(1 + k*r))/(2*k)
    return y, dy

def extend_Met1(En, Mas, rD, sD, dsD, uD, duD, Ext,
                Np=1000):
    """
    Extends the radial profile of a field using asymptotic approximations.

    Parameters:
        En (float): Energy parameter.
        Mas (float): Mass parameter.
        rD (np.ndarray): Array of radial points.
        sD (np.ndarray): Scalar field values at rD.
        dsD (np.ndarray): Derivative of scalar field at rD.
        uD (np.ndarray): Function U values at rD.
        duD (np.ndarray): Derivative of U at rD.
        Ext (float): Extension length for the radial coordinate.
        Np (int, optional): Number of new points for extension (default: 1000).

    Returns:
        tuple: Extended arrays (rDnew, sDnew, dsDnew, uDnew, duDnew).
    """

    # Validate input
    if len(rD) == 0 or np.isnan(rD[-1]):
        raise ValueError("rD must be a non-empty array with valid numerical values.")

    # Generate extended radial points
    rad = np.linspace(rD[-1], rD[-1] + Ext, Np)

    # Compute asymptotic parameters
    Ap, k = sParam(rD, sD)

    # Compute asymptotic profiles
    sExt, dsExt = Asy_sProf(rad, Ap, k)
    uExt, duExt = Asy_uProf(rad, En, Mas)

    # Concatenate with existing data
    rDnew = np.concatenate((rD[:-1], rad))
    sDnew = np.concatenate((sD[:-1], sExt))
    dsDnew = np.concatenate((dsD[:-1], dsExt))
    uDnew = np.concatenate((uD[:-1], uExt))
    duDnew = np.concatenate((duD[:-1], duExt))

    return rDnew, sDnew, dsDnew, uDnew, duDnew

def extend_Met2(En, Mas, rD, sD, dsD, uD, duD, Ext,
                Np=1000, info=False):
    """
    Extends the radial profile of a field using asymptotic approximations.

    Parameters:
        En (float): Energy parameter.
        Mas (float): Mass parameter.
        rD (np.ndarray): Array of radial points.
        sD (np.ndarray): Scalar field values at rD.
        dsD (np.ndarray): Derivative of scalar field at rD.
        uD (np.ndarray): Function U values at rD.
        duD (np.ndarray): Derivative of U at rD.
        Ext (float): Extension length for the radial coordinate.
        Np (int, optional): Number of new points for extension (default: 1000).
        info: Print the energy compute by two via

    Returns:
        tuple: Extended arrays (rDnew, sDnew, dsDnew, uDnew, duDnew).
    """

    # Validate input
    if len(rD) == 0 or np.isnan(rD[-1]):
        raise ValueError("rD must be a non-empty array with valid numerical values.")

    # Generate extended radial points
    rad = np.linspace(rD[-1], rD[-1] + Ext, Np)

    # Compute asymptotic parameters
    Ap, k = sParam(rD, sD, B=Mas)
    if info:
        print('En = ', En, ' from k => En =', -k**2)

    # Compute asymptotic profiles
    sExt, dsExt = Asy_sProf_v2(rad, Ap, En, Mas)
    uExt, duExt = Asy_uProf(rad, En, Mas)

    # Concatenate with existing data
    rDnew = np.concatenate((rD[:-1], rad))
    sDnew = np.concatenate((sD[:-1], sExt))
    dsDnew = np.concatenate((dsD[:-1], dsExt))
    uDnew = np.concatenate((uD[:-1], uExt))
    duDnew = np.concatenate((duD[:-1], duExt))

    return rDnew, sDnew, dsDnew, uDnew, duDnew

__all__ = ['mainExt', 'extend', 'sParam', 'Asy_uProf', 'Asy_sProf', 'Asy_sProf_v2', 'extend_Met1', 'extend_Met2']
