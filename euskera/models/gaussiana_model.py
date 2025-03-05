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


########### Class for the gaussiana model
#############################################################################

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
        dat_conf_Parameters = zip(*(self.parameters_gaussiana[i] for i in self.parameters_gaussiana.keys()))
        
        num_threads = parameters_simulation.get("num_threads", 1)
        ne.set_num_threads(num_threads)  # Set the number of threads for parallel execution
        set_num_threads(num_threads)  # Set the number of threads for parallel execution of numba
        
        if psi is not None and xarray is not None:
            for posGauss, amp, sig in dat_conf_Parameters:
                psi = initGaussianInt(psi, xarray, yarray, zarray, posGauss, sig, amp=1.0 if amp is None else amp)
            # Compute the density profile for every component
            rho_i = ne.evaluate("real(abs(psi)**2)")
            return [psi, rho_i]
        else:
            # Extract simulation parameters from dictionary
            resol = parameters_simulation.get("resol", 128)
            gridlength = parameters_simulation.get("gridlength", 1.0)
            num_threads = parameters_simulation.get("num_threads", 1)
            
            # Generate spatial grids and distance array
            [xarray, yarray, zarray], distarray = gd.RealGrid(gridlength=gridlength, resol=resol)
            
            # Create aligned zero-filled arrays for Psi components
            # Use pyFFTW if available, otherwise fallback to NumPy
            if pyfftwOpt:
                psi = pyfftw.zeros_aligned((field_components, resol, resol, resol), dtype='complex128')
            else:
                psi = np.zeros((field_components, resol, resol, resol), dtype='complex128')
                
            for posGauss, amp, sig in dat_conf_Parameters:
                psi = initGaussianInt(psi, xarray, yarray, zarray, posGauss, sig, amp=1.0 if amp is None else amp)
            
            # Compute the density profile for every component
            rho_i = ne.evaluate("real(abs(psi)**2)")
            return [xarray, yarray, zarray, distarray], [psi, rho_i]
        

########### Extra functions (Outside of the clase)
#############################################################################
@njit(fastmath=True)
def gaussian(dx, dy, dz, sigx, sigy, sigz, amp):
    """
    Computes the Gaussian function value for a given coordinate.
    Optimized with Numba for performance.
    """
    fg = amp * np.exp(-(dx**2 / (2 * sigx**2) +
                          dy**2 / (2 * sigy**2) +
                          dz**2 / (2 * sigz**2)))
    return fg

@njit(parallel=True, fastmath=True)
def initGaussianInt(funct, xarray, yarray, zarray, posGauss, sig, amp=1.0):
    """
    Initializes a 3D Gaussian intensity function over a given grid.
    Optimized with Numba for parallel execution.
    """
    sigx, sigy, sigz = sig  # Unpack sigma values
    amp = float(amp)

    dx = xarray[:, 0, 0] - posGauss[0]
    dy = yarray[0, :, 0] - posGauss[1]
    dz = zarray[0, 0, :] - posGauss[2]

    for sol in prange(funct.shape[0]):  # Use parallel execution
        for i in prange(funct.shape[1]):
            for j in range(funct.shape[2]):
                for k in range(funct.shape[3]):
                    funct[sol, i, j, k] += gaussian(dx[i], dy[j], dz[k], sigx, sigy, sigz, amp)

    return funct
    
#@njit(parallel=True)
#def initGaussianInt(funct, xarray, yarray, zarray, posGauss, sig, amp=1.0):
#    """
#    Initializes a 3D Gaussian intensity function over a given grid.
#    Optimized with Numba for parallel execution.
#    """
#    sigx, sigy, sigz = sig  # Unpack sigma values
#    amp = float(amp)
#    
#    for sol in prange(funct.shape[0]):  # Use parallel execution
#        for i in range(funct.shape[1]):
#            for j in range(funct.shape[2]):
#                for k in range(funct.shape[3]):
#                    dx = xarray[i, 0, 0] - posGauss[0]
#                    dy = yarray[0, j, 0] - posGauss[1]
#                    dz = zarray[0, 0, k] - posGauss[2]
#                
#                    funct[sol, i, j, k] += gaussian(dx, dy, dz, sigx, sigy, sigz, amp)
#    return funct