# euskera v1.0
# soliton_model file

import sys
import os
import numpy as np
import numexpr as ne

from numba import njit, prange  

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
class Soli_Modelo():
    """ 
    Class used for modeling solitons
    """
    
    # Class attribute
    pyfftwOpt = pyfftwOpt  # Store FFTW availability as a class attribute
    
    def __init__(self, parameters_mod):
        self.parameters_mod = parameters_mod   
    
    def SoliIncProf(self, psi, coord, field_components, parameters_simulation):
        """
        Computes the initial soliton profiles and updates the wavefunction psi.

        Args:
            psi (np.ndarray): Initial wavefunction array.
            coord (list): List containing x, y, and z coordinate arrays.
            field_components (int): Number of field components.
            parameters_simulation (list): Simulation parameters.

        Returns:
            np.ndarray: Updated wavefunction array with soliton profiles.
        """
        name_param_by_configuration = {
            "positions": 3, "velocities": 3, "phases": 1,  "alphas": 1, "dr": 1
        }
        name_param_by_components = {
            "profiles": field_components, "betas": field_components
            }
        
        modelo_parameters = self.parameters_mod
        ###################################################################################################
        # Validate configuration parameters
        for name, expected_quantity in name_param_by_configuration.items():
            parameter_list = modelo_parameters.get(name, [])
            for element in parameter_list:  
                actual_quantity = len(element)
                if actual_quantity != expected_quantity:
                    raise ValueError(
                        f"The parameter '{name}' must have {expected_quantity} components, but got {actual_quantity}.")

        # Validate component parameters
        for name, expected_quantity in name_param_by_components.items():
            parameter_list = modelo_parameters.get(name, [])
            for element in parameter_list:
                actual_quantity = len(element)
                if actual_quantity != expected_quantity:
                    raise ValueError(f"The parameter '{name}' must have {expected_quantity} elements, but got {actual_quantity}.")
        ###################################################################################################
        
        # Extracting simulation parameters
        Boverlap, Plim, rmax, resol, lambda_value, t0 = parameters_simulation
        
        if lambda_value != 0:
            print("WARNING: The alpha values are fixed to 1 because the scaling property is not true.")
            modelo_parameters["phases"] = [[1.]]*len(modelo_parameters["phases"])
        ###################################################################################################
        
        # Handle soliton overlap check if enabled
        if Boverlap:
            # Check input validity
            if Plim >= rmax:
                raise ValueError(f"Plim ({Plim}) must be less than rmax ({rmax}).")
    
            warn = to.overlap(modelo_parameters["positions"], rmax=rmax, warn=0)
            if warn:
                print("WARNING: Significant overlap detected between solitons in initial conditions \n")
        ###################################################################################################
                
        # Zip the configuration parameters
        conf_parameters_name = list(name_param_by_configuration.keys())
        dat_conf_Parameters = zip(*(modelo_parameters[i] for i in conf_parameters_name))
        
        # Zip the component parameters with configuration parameters
        dat_glob_Parameters = zip(
            modelo_parameters["profiles"],
            modelo_parameters["betas"],
            dat_conf_Parameters
        )
        
        # Create aligned zero-filled arrays for the temporal funct components
        # Use pyFFTW if available, otherwise fallback to NumPy
        if self.pyfftwOpt:
            funct = pyfftw.zeros_aligned((resol, resol, resol), dtype='complex128')
        else:
            funct = np.zeros((resol, resol, resol), dtype='complex128')
        
        # Generate the initial soliton profiles
        xarray, yarray, zarray = coord
        for comp_profiles, beta_comp_values, conf_param in dat_glob_Parameters:
            positionCen, comp_veloc, phase, alpha, dr = conf_param
            for i in range(field_components):  # running by the components
                psi_i = psi[i]
                fInt = comp_profiles[i]
                
                param = [beta_comp_values[i], phase[0], positionCen, alpha[0]]
                psi[i] = InitialSolitonsProf(
                    psi_i, funct, fInt, [xarray, yarray, zarray],
                    comp_veloc, param, Plim=Plim, t0=t0, delta_x=dr[0],
                    initsoliton=initsolitonInt
                )          
        return psi
    
    def PsiInic(self, field_components, parameters_simulation):
        """
        Initializes a soliton field in a 3D grid.

        Args:
            field_components (int): Number of field components.
            parameters_simulation (dict): Simulation parameters.

        Returns:
            list: Spatial grids, distance array, wavefunction, and density profile.
        """
        
        # Extract simulation parameters from dictionary
        Boverlap = parameters_simulation.get("Boverlap", False)
        Plim = parameters_simulation.get("Plim", 1.0)
        rmax = parameters_simulation.get("rmax", 1.0)
        resol = parameters_simulation.get("resol", 128)
        lambda_value = parameters_simulation.get("lambda_value", 1.0)
        t0 = parameters_simulation.get("t0", 0.0)
        gridlength = parameters_simulation.get("gridlength", 1.0)
        num_threads = parameters_simulation.get("num_threads", 1)

        # Generate spatial grids and distance array
        [xarray, yarray, zarray], distarray = gd.RealGrid(gridlength=gridlength, resol=resol)
    
        # Create aligned zero-filled arrays for Psi components
        # Use pyFFTW if available, otherwise fallback to NumPy
        if Soli_Modelo.pyfftwOpt:
            psi = pyfftw.zeros_aligned((field_components, resol, resol, resol), dtype='complex128')
        else:
            psi = np.zeros((field_components, resol, resol, resol), dtype='complex128')

        # Compute the initial soliton profile
        parameters_simulation2 = [Boverlap, Plim, rmax, resol, lambda_value, t0]
        coord = [xarray, yarray, zarray]
        psi = self.SoliIncProf(psi, coord, field_components, parameters_simulation2)

        # Set the number of threads for parallel execution
        ne.set_num_threads(num_threads)

        # Compute the density profile for every component
        rho_i = ne.evaluate("real(abs(psi)**2)")
    
        return [xarray, yarray, zarray, distarray], [psi, rho_i]
 

########### Extra functions (Outside of the clase)
#############################################################################
def InitialSolitonsProf(psi, funct, fInt, grid, velocity, param, initsoliton, Plim=5.6,
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
    