"""SOLITON MODEL FOR EUSKERA SIMULATIONS"""

import numpy as np
import numexpr as ne

from numba import njit, prange, set_num_threads

try:
    import pyfftw
    pyfftwOpt = True
except ImportError:
    pyfftwOpt = False  


from euskera.core import grids as gd
from euskera.tools import tools as to
###################################################################################################


# ============================================================
# ====================== SOLITON MODEL =======================
# ============================================================
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

        shape = (resol, resol, resol)
        if field_components: 
            shape = (field_components,) + shape

        out = pyfftw.zeros_aligned(shape, dtype='complex128') if self.pyfftwOpt else np.zeros(shape, dtype='complex128')
        return out

    def add_soliton(self, psi, coord, field_components, parameters_simulation):
        """Computes the initial soliton profiles and updates psi."""
        self.validate_parameters(field_components)

        Boverlap, Plim, rmax, resol, lambda_value, num_threads = parameters_simulation
        set_num_threads(num_threads)
        
        if lambda_value != 0:
            print("WARNING: Alpha values are fixed to 1 due to scaling constraints.")
            self.parameters_mod["alphas"] = [[1.0]] * len(self.parameters_mod["alphas"])
        
        if Boverlap and Plim >= rmax:
            raise ValueError(f"Plim ({Plim}) must be less than rmax ({rmax}).")
        
        if Boverlap and to.overlap(self.parameters_mod["positions"], rmax=rmax, warn=0):
            print("WARNING: Significant soliton overlap detected.")
        
        funct = self.initialize_wavefunction(resol)

        # Evalueate the Psi function over the grid
        # Zip the configuration parameters
        conf_names = ["positions", "velocities", "phases", "alphas", "dr"]
        conf_data = zip(*(self.parameters_mod[name] for name in conf_names))

        # Zip the component parameters with configuration parameters
        glob_data = zip(self.parameters_mod["profiles"], self.parameters_mod["betas"], conf_data)

        for comp_profiles, beta_values, conf in glob_data:
            position, velocity, phase, alpha, dr = conf

            # Loop over vector components
            for i in range(field_components):
                param = (float(beta_values[i]), float(phase[0]), position, float(alpha[0]))
                psi[i] = build_soliton(funct, psi[i], comp_profiles[i], coord, velocity, param, Plim=Plim, delta_x=float(dr[0]))
        return psi
    

    def apply(self, field_components, parameters_simulation, grid_data=None, psi=None):
        """Initializes an soliton-boson field in a 3D grid"""

        if (psi is None) != (grid_data is None): raise ValueError("psi and grid arrays must be both provided or both None")

        build_grid = psi is None
        
        # Grid generation
        if build_grid:
            # Extract simulation parameters from dictionary
            resol = parameters_simulation.get("resol", 128)
            gridlength = parameters_simulation.get("gridlength", 1.0)

            # Generate spatial grids and distance array
            [xarray, yarray, zarray], distarray = gd.RealGrid(gridlength=gridlength, resol=resol)
            grid_data = [xarray, yarray, zarray, distarray]

            # Allocate wavefunction
            psi = self.initialize_wavefunction(parameters_simulation["resol"], field_components)
        else:
            xarray, yarray, zarray, distarray = grid_data
        
        # Apply Soliton
        # Puting the soliton profiles on Psi
        coord = grid_data[:3]
        sim_params = [parameters_simulation["Boverlap"], parameters_simulation["Plim"], parameters_simulation["rmax"],
                      parameters_simulation["resol"], parameters_simulation["lambda_value"], parameters_simulation["num_threads"]
                    ]
        
        psi = self.add_soliton(psi, coord, field_components, sim_params)
        
        # Set the number of threads for parallel execution
        ne.set_num_threads(parameters_simulation["num_threads"])
        
        # Compute the density profile for every component   
        rho_i = ne.evaluate("real(psi * conj(psi))") # ("real(abs(psi)**2)")
        
        return grid_data, psi, rho_i
 

########### Auxiliary Functions Outside of the class)
#############################################################################
             
def build_soliton(funct, psi, radial_profile, 
                  grid, velocity, param, 
                  Plim=5.6, delta_x=1e-5, t0=0):
    """
    Compute the initial soliton profile:
    
    Parameters:
    - psi: Configuration Wavefunction.
    - funct: wave function
    - fInt: Function that computes the soliton profile.
    - grid: Tuple containing (xarray, yarray, zarray).
    - velocity: Tuple (velx, vely, velz).
    - param: Tuple (beta, phase, positionCen, alpha).
    - Plim: Utilises soliton profile array out to dimensionless radius (default: 5.6).
    
    Returns:
    - psi: Wavefunction.
    - funct: Soliton profile before multiuply by the factor: exp(i(alpha beta t +vect{v}*vect{x}-1/2|v|^2 t))
    """
    
    beta, phase, positionCen, alpha = param 
    xarray, yarray, zarray = grid
    velx, vely, velz = velocity

    # Radial
    funct = initsoliton_kernel(funct, xarray, yarray, zarray, positionCen, radial_profile, Plim=Plim, alpha=alpha, delta_x=delta_x)
    
    ####### Impart velocity to solitons in Galilean invariant way
    funct = ne.evaluate("exp(1j*(alpha*beta*t0 + velx*xarray + vely*yarray + velz*zarray - 0.5*(velx*velx + vely*vely + velz*velz)*t0 + phase)) * funct")
    psi = ne.evaluate("psi + funct")
    
    return psi

@njit(parallel=True)
def initsoliton_kernel(funct, xarray, yarray, zarray, position,
                   radial_profile, Plim=5.6, alpha=1., delta_x=1e-5):
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
    rmax_sq = Plim * Plim
    inv_delta = 1.0 / delta_x
    max_index = len(radial_profile) - 1

    nx, ny, nz = funct.shape
    
    for i in prange(nx):  # Use parallel execution
        dx = xarray[i, 0, 0] - position[0]
        for j in range(ny):
            dy = yarray[0, j, 0] - position[1]
            for k in range(nz):
                dz = zarray[0, 0, k] - position[2]
                
                dist_sq = dx * dx + dy * dy + dz * dz  # Squared distance
                scaled_dist_sq = alpha * dist_sq if alpha is not None else dist_sq

                if scaled_dist_sq <= rmax_sq:
                    r = np.sqrt(scaled_dist_sq)
                    index2 = int(r * inv_delta)

                    if index2 > max_index:
                        index2 = max_index

                    funct[i, j, k] = (
                        alpha * radial_profile[index2]
                        if alpha is not None
                        else radial_profile[index2]
                    )
                else:
                    funct[i, j, k] = 0.0 

    return funct