# euskera v1.0
# soliton_model file

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


# Get the parent directory dynamically
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parent_dir)

import main.grids as gd
import tools.tools as to
###################################################################################################


########### Class for the solitonic model
#############################################################################
class Soli_Model():
    """ 
    Class used for modeling solitons
    """
    
    # Class attribute
    pyfftwOpt = pyfftwOpt  # Store FFTW availability as a class attribute
    
    def __init__(self, parameters_mod):
        self.parameters_mod = parameters_mod
    
    def validate_parameters(self, field_components):
        """Validates the input parameters for configuration and component correctness."""
        param_config = {"positions": 3, "velocities": 3, "phases": 1, "alphas": 1, "dr": 1}
        param_components = {"profiles": field_components, "betas": field_components}
        
        # Validate configuration parameters
        for name, expected in param_config.items():
            for element in self.parameters_mod.get(name, []):
                if len(element) != expected:
                    raise ValueError(f"Parameter '{name}' must have {expected} elements, but got {len(element)}.")
        
        # Validate component parameters
        for name, expected in param_components.items():
            for element in self.parameters_mod.get(name, []):
                if len(element) != expected:
                    raise ValueError(f"Parameter '{name}' must have {expected} elements, but got {len(element)}.")
    
    def initialize_wavefunction(self, resol, field_components=None):
        """Allocates memory for the wavefunction psi."""
        # Create aligned zero-filled arrays for the temporal funct components
        # Use pyFFTW if available, otherwise fallback to NumPy
        
        if field_components:
            out = pyfftw.zeros_aligned((field_components, resol, resol, resol), dtype='complex128') if self.pyfftwOpt \
                else np.zeros((field_components, resol, resol, resol), dtype='complex128')
        else:
            out = pyfftw.zeros_aligned((resol, resol, resol), dtype='complex128') if self.pyfftwOpt \
                else np.zeros((resol, resol, resol), dtype='complex128')
        return out

    def SoliIncProf(self, psi, coord, field_components, parameters_simulation):
        """Computes the initial soliton profiles and updates psi."""
        self.validate_parameters(field_components)
        Boverlap, Plim, rmax, resol, lambda_value, t0, num_threads = parameters_simulation
        
        if lambda_value != 0:
            print("WARNING: Alpha values are fixed to 1 due to scaling constraints.")
            self.parameters_mod["alphas"] = [[1.]] * len(self.parameters_mod["alphas"])
        
        if Boverlap and Plim >= rmax:
            raise ValueError(f"Plim ({Plim}) must be less than rmax ({rmax}).")
        
        if Boverlap and to.overlap(self.parameters_mod["positions"], rmax=rmax, warn=0):
            print("WARNING: Significant overlap detected between solitons.")
        
        funct = self.initialize_wavefunction(resol)
        psi = psievaluation(funct, psi, self.parameters_mod,
                            initsolitonInt, coord, field_components, Plim, t0, num_threads)
        return psi
    
    
    def PsiInic(self, field_components, parameters_simulation):
        """Initializes a soliton field in a 3D grid."""
        coord, distarray = gd.RealGrid(gridlength=parameters_simulation["gridlength"], resol=parameters_simulation["resol"])
        psi = self.initialize_wavefunction(parameters_simulation["resol"], field_components)
        
        # Generate the initial soliton profiles
        parameters_simulation2 = [parameters_simulation["Boverlap"], parameters_simulation["Plim"], parameters_simulation["rmax"],
                                  parameters_simulation["resol"], parameters_simulation["lambda_value"], parameters_simulation["t0"],
                                  parameters_simulation["num_threads"]]
        psi = self.SoliIncProf(psi, coord, field_components, parameters_simulation2)
        
        # Set the number of threads for parallel execution
        ne.set_num_threads(parameters_simulation["num_threads"])
        # Compute the density profile for every component   
        rho_i = ne.evaluate("real(abs(psi)**2)")
        return coord + [distarray], [psi, rho_i]
 

########### Auxiliary Functions Outside of the class)
#############################################################################

def psievaluation(funct, psi, params, initsolitonInt, coord, field_components, Plim, t0, num_threads):
    
    # Zip the configuration parameters
    conf_parameters_name = ["positions", "velocities", "phases", "alphas", "dr"]
    dat_conf_Parameters = zip(*(params[i] for i in conf_parameters_name))
        
    # Zip the component parameters with configuration parameters
    dat_glob_Parameters = zip(params["profiles"], params["betas"], dat_conf_Parameters)
        
    set_num_threads(num_threads)
    for comp_profiles, beta_comp_values, conf_param in dat_glob_Parameters:
        positionCen, comp_veloc, phase, alpha, dr = conf_param
        for i in range(field_components): # running by the components
            param = [beta_comp_values[i], phase[0], positionCen, alpha[0]]
            psi[i] = InitialSolitonsProf(funct, psi[i], comp_profiles[i],
                                         coord, comp_veloc, param, Plim=Plim, t0=t0, delta_x=dr[0], initsoliton=initsolitonInt)
    return psi
             
def InitialSolitonsProf(funct, psi, fInt, grid, velocity, param, initsoliton, Plim=5.6,
                        t0=0, delta_x=0.00001):
    """
    Compute the initial soliton profile:
    
    \psi_i(\vec{x}, t) = scal fInt(\sqrt{scal}|\vec{x}-\vec{v}t|) 
                         \exp(i(- scal \beta t + \vec{v} \cdot \vec{x} - |\vec{v}|^2 t/2))
    
    Parameters:
    - psi: Configuration Wavefunction.
    - funct: wave function
    - fInt: Function that computes the soliton profile.
    - grid: Tuple containing (xarray, yarray, zarray).
    - velocity: Tuple (velx, vely, velz).
    - param: Tuple (beta, phase, positionCen, alpha).
    - Plim: Utilises soliton profile array out to dimensionless radius (default: 5.6).
    - t0: Initial time (default 0).
    
    Returns:
    - psi: Wavefunction.
    - funct: Soliton profile before multiuply by the factor: exp(i(alpha beta t +\vect{v}*\vect{x}-1/2|v|^2 t))
    """
    
    beta, phase, positionCen, alpha = param 
    xarray, yarray, zarray = grid
    velx, vely, velz = velocity

    # Compute initial soliton shape
    if t0 == 0:
        funct = initsoliton(funct, grid[0], grid[1], grid[2],
                            positionCen, fInt, Plim=Plim, alpha=alpha, delta_x=delta_x)
    else:
        pass # poner el caso con dx - v*t

    ####### Impart velocity to solitons in Galilean invariant way
    funct = ne.evaluate("exp(1j*(alpha*beta*t0 + velx*xarray + vely*yarray + velz*zarray - 0.5*(velx*velx + vely*vely + velz*velz)*t0 + phase)) * funct")
    psi = ne.evaluate("psi + funct")
    return psi


@njit(parallel=True)
def initsolitonInt(funct, xarray, yarray, zarray, position,
                   fInt, Plim=5.6, alpha=1., delta_x=0.00001):
    """
    Initializes a soliton profile in a 3D grid using Numba optimization.
    
    Parameters:
    - funct: 3D NumPy array to store the soliton function values.
    - xarray, yarray, zarray: 3D NumPy arrays representing spatial coordinates.
    - position: NumPy array (x, y, z) indicating the soliton center.
    - fInt: Tuple (dr, sigV) for interpolation.
    - rmax: Maximum radius for the soliton (default 5.6).
    - alpha: Scaling factor (default None).

    Returns:
    - Updated `funct` array with soliton values.
    """

    rmax_sq = Plim ** 2  # Precompute squared max radius
    
    for i in prange(funct.shape[0]):  # Use parallel execution
        for j in range(funct.shape[1]):
            for k in range(funct.shape[2]):
                dx = xarray[i, 0, 0] - position[0]
                dy = yarray[0, j, 0] - position[1]
                dz = zarray[0, 0, k] - position[2]
                
                dist_sq = dx * dx + dy * dy + dz * dz  # Squared distance
                
                if alpha is not None:
                    scaled_dist_sq = alpha * dist_sq
                    if scaled_dist_sq <= rmax_sq:
                        index2 = int(np.sqrt(scaled_dist_sq) / delta_x)
                        funct[i, j, k] = alpha * fInt[index2]
                    else:
                        funct[i, j, k] = 0
                else:
                    if dist_sq <= rmax_sq:
                        index2 = int(np.sqrt(dist_sq) / delta_x)
                        funct[i, j, k] = fInt[index2]
                    else:
                        funct[i, j, k] = 0

    return funct

def initsoliton(funct, xarray, yarray, zarray, position, fInt, Plim=5.6, delta_x=None, alpha=None):
    """
    Initializes a soliton profile in a 3D grid. Without numba optimization
    
    Note that we compute the distance of every gridpoint from the centre of the soliton, 
    not to calculate the distance of the soliton from the centre of the grid

    Parameters:
    - funct: 3D NumPy array to store the soliton function values.
    - xarray, yarray, zarray: 3D NumPy arrays representing spatial coordinates.
    - position: Tuple or array (x, y, z) indicating the soliton center.
    - fInt: Function that computes the soliton profile.
    - Plim: Utilises soliton profile array out to dimensionless radius (default: 5.6).
    - alpha: Scaling factor (default None).

    Returns:
    - Updated `funct` array with soliton values.
    """
    
    if not isinstance(position, (list, tuple, np.ndarray)) or len(position) != 3:
        raise ValueError("position must be a list, tuple, or NumPy array of length 3.")
    
    for index in np.ndindex(funct.shape):
        dx = xarray[index[0], 0, 0] - position[0]
        dy = yarray[0, index[1], 0] - position[1]
        dz = zarray[0, 0, index[2]] - position[2]
        
        dist_sq = dx**2 + dy**2 + dz**2  # Squared distance
        if alpha:
            scaled_dist_sq = alpha * dist_sq
            if scaled_dist_sq <= Plim**2:  # Compare squared values to avoid sqrt
                funct[index] = alpha * fInt(np.sqrt(scaled_dist_sq))
            else:
                funct[index] = 0
        else:
            if dist_sq <= Plim**2:
                funct[index] = fInt(np.sqrt(dist_sq))
            else:
                funct[index] = 0
        
    return funct
    