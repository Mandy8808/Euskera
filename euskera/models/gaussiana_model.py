# euskera v1.0
# gaussiana_model file

import sys
import os
import numpy as np
import numexpr as ne 

from numba import njit, prange, set_num_threads

try:
    import pyfftw
    pyfftwOpt = True
except ImportError:
    pyfftwOpt = False

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ""))
sys.path.append(parent_dir)

import main.grids as gd

###################################################################################################

# ============================================================
# ==================== GAUSSIANA MODEL =======================
# ============================================================

class Gaussiana_Model():
    """ 
    Class for modeling Gaussian perturbations in field components
    """
    
    # Class attribute
    pyfftwOpt = pyfftwOpt  # Store FFTW availability as a class attribute
    
    def __init__(self, parameters_mod):
        self.parameters_gaussiana = parameters_mod
        
    def GaussSolProf(self, field_components, parameters_simulation,
                    psi=None, xarray=None, yarray=None, zarray=None):
        """
        Applies Gaussian perturbations and optionally Soliton Profile
        """
        if (psi is None) != (xarray is None):
            raise ValueError("psi and grid arrays must be both provided or both None")
        
        num_threads = parameters_simulation.get("num_threads", 1)
        ne.set_num_threads(num_threads)  # Set the number of threads for parallel execution
        set_num_threads(num_threads)  # Set the number of threads for parallel execution of numba

        dat_conf_Parameters = zip(*(self.parameters_gaussiana[i] for i in self.parameters_gaussiana.keys()))

        build_grid = psi is None
        if build_grid:
            # Extract simulation parameters from dictionary
            resol = parameters_simulation.get("resol", 128)
            gridlength = parameters_simulation.get("gridlength", 1.0)

            # Generate spatial grids and distance array
            [xarray, yarray, zarray], distarray = gd.RealGrid(gridlength=gridlength, resol=resol)

            # Create aligned zero-filled arrays for Psi components
            # Use pyFFTW if available, otherwise fallback to NumPy   
            shape = (field_components, resol, resol, resol)
            psi = (pyfftw.zeros_aligned(shape, dtype='complex128')
                   if self.pyfftwOpt
                   else np.zeros(shape, dtype='complex128')
            )
        
        # Puting the Gaussinan profiles on Psi
        for posGauss, amp, sig in dat_conf_Parameters:
            amp = 1.0 if amp is None else amp
            psi = initGaussianInt(psi, xarray, yarray, zarray, posGauss, sig, amp=amp)

        # Compute the density profile for every component
        rho_i = ne.evaluate("real(abs(psi)**2)")

        if build_grid:
            return [xarray, yarray, zarray, distarray], [psi, rho_i]
        else:
            return [psi, rho_i]
        

########### Extra functions (Outside of the clase)
#############################################################################
@njit(fastmath=True)
def build_1d_gaussian(coord, center, sigma):
    """Computes the Gaussian function"""
    out = np.empty_like(coord)
    inv = 1.0 / (2.0 * sigma * sigma)
    for i in range(coord.shape[0]):
        dx = coord[i] - center
        out[i] = np.exp(-dx * dx * inv)
    return out

@njit(parallel=True, fastmath=True)
def initGaussianInt(funct, xarray, yarray, zarray, posGauss, sig, amp=1.0):
    """
    Initializes a 3D Gaussian intensity function over a given grid.
    Optimized with Numba for parallel execution.
    """
    sigx, sigy, sigz = sig  # Unpack sigma values

    # Extraer ejes 1D
    x = xarray[:, 0, 0]
    y = yarray[0, :, 0]
    z = zarray[0, 0, :]

    gx = build_1d_gaussian(x, posGauss[0], sigx)
    gy = build_1d_gaussian(y, posGauss[1], sigy)
    gz = build_1d_gaussian(z, posGauss[2], sigz)

    amp = float(amp)

    for comp in prange(funct.shape[0]):  # Use parallel execution
        for i in range(funct.shape[1]):
            for j in range(funct.shape[2]):
                for k in range(funct.shape[3]):
                    funct[comp, i, j, k] += amp * gx[i] * gy[j] * gz[k]

    return funct