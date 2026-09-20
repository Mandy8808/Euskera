"""GAUSSIANA MODEL FOR EUSKERA SIMULATIONS"""

import numpy as np
import numexpr as ne

from numba import njit, prange, set_num_threads

try:
    import pyfftw
    pyfftwOpt = True
except ImportError:
    pyfftwOpt = False

from euskera.core import grids as gd

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
    
    def initialize_wavefunction(self, resol, field_components=None):
        """Allocates memory for the wavefunction psi."""
        # Create aligned zero-filled arrays for the temporal funct components
        # Use pyFFTW if available, otherwise fallback to NumPy

        shape = (resol, resol, resol)
        if field_components: 
            shape = (field_components,) + shape

        out = pyfftw.zeros_aligned(shape, dtype='complex128') if self.pyfftwOpt else np.zeros(shape, dtype='complex128')
        return out
        
    def apply(self, field_components, parameters_simulation, grid_data=None, psi=None):
        """
        Applies Gaussian perturbations and optionally Soliton Profile
        """
        if (psi is None) != (grid_data is None): raise ValueError("psi and grid arrays must be both provided or both None")
        
        num_threads = parameters_simulation.get("num_threads", 1)
        ne.set_num_threads(num_threads)  # Set the number of threads for parallel execution
        set_num_threads(num_threads)  # Set the number of threads for parallel execution of numba

        build_grid = psi is None
        if build_grid:
            # Extract simulation parameters from dictionary
            resol = parameters_simulation.get("resol", 128)
            gridlength = parameters_simulation.get("gridlength", 1.0)

            # Generate spatial grids and distance array
            [xarray, yarray, zarray], distarray = gd.RealGrid(gridlength=gridlength, resol=resol)
            grid_data = [xarray, yarray, zarray, distarray]

            # Allocate wavefunction
            psi = self.initialize_wavefunction(resol, field_components)
        else:
            xarray, yarray, zarray, distarray = grid_data
        
        # Apply Gaussians
        # Puting the Gaussinan profiles on Psi
        dat_conf_Parameters = zip(*(self.parameters_gaussiana[i] for i in self.parameters_gaussiana.keys()))
        
        for posGauss, amp, sig in dat_conf_Parameters:
            amp = 1.0 if amp is None else amp
            posGauss = np.asarray(posGauss, dtype=np.float64)
            sig = np.asarray(sig, dtype=np.float64)
            psi = add_gaussian(psi, xarray, yarray, zarray, posGauss, sig, amp=amp)

        # Compute the density profile for every component
        rho_i = ne.evaluate("real(psi * conj(psi))") # ne.evaluate("real(abs(psi)**2)")

        return grid_data, psi, rho_i
        

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
def add_gaussian(funct, xarray, yarray, zarray, posGauss, sig, amp=1.0):
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