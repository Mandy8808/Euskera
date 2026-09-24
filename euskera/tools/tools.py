# euskera v1.0
# tools file

import numpy as np
import os

from scipy.integrate import quad
from scipy.interpolate import interp1d

#############################################################################

########### Functions to updated default parameters
#############################################################################
def update_simulation_parameters(simulation_parameters_update, simulation_parameters):
    """
    Updates simulation parameters with validation.

    Args:
        simulation_parameters_update (dict): Parameters to update.
        simulation_parameters (dict): Existing parameters dictionary.

    Returns:
        dict: Updated simulation parameters.

    Raises:
        ValueError: If an update has an incorrect type.
    """
    
    if not isinstance(simulation_parameters_update, dict):
        raise ValueError("The simulation_parameters_update must be of type dict.")
    
    if not isinstance(simulation_parameters, dict):
        raise ValueError("The simulation_parameters must be of type dict.")

    simulation_type = {
        "lambda_value": (int, float),
        "num_threads": int,
        "gridlength": (int, float),
        "resol": int,
        "step_factor": (int, float),
        "t0": (int, float),
        "tmax": (int, float),
        "rmax": (int, float),
        "Plim": (int, float),
        "cmass": (int, float),
        "plott0": bool,
        "Boverlap": bool,
        "info": bool,
        "save_number": int,
        "cleanup_policy": str,
        "consolidation_batch_size": int,
        "format": str,
        "address": str,
        "copy_profiles": bool,
        "data_save": dict,
        "rules": (list, tuple, type(None)),
        "Numb_Part": bool,
        "Energ": bool,
        "Pi": bool,
        "Ji": bool,
        "methodEnerg": int,
        "Frequency": bool
    }

    for name, element_update in simulation_parameters_update.items():
        expected_type = simulation_type.get(name)

        if expected_type is None:
            raise ValueError(f"Unknown parameter: '{name}' is not a valid simulation parameter.")

        if not isinstance(element_update, expected_type):
            raise ValueError(f"The parameter '{name}' must be of type {expected_type}.")
        if name == "Ji" and element_update:
            raise ValueError(
                "Ji diagnostics are not implemented; set Ji=False."
            )

        simulation_parameters[name] = element_update

    return simulation_parameters


########### Extract time step information
#############################################################################
def dtime(tmax, gridlength, resol, step_factor, save_number):
    """
    Computes the size of the timestep (∆t). (CAN BE INCREASED WITH step_factor)
    
    Comments:
        The timestep ∆t, by default is take so that fluid travelling at this maximum velocity traverses one grid space, ∆x, per timestep:
            ∆t = (∆x)^2/π

    Parameters:
        tmax (float): Duration time in code units.
        gridlength (float): The physical length of the grid in one dimension.
        resol (int): The number of points in each dimension (resolution).
        step_factor (float): Factor to increase timestep if velocities are low.
        save_number (int): Number of frames to save.

    Returns:
        ht (float): Time step per iteration.
        its_per_save (int): Iterations per save step.
    """
    
    values = {
        "tmax": tmax,
        "gridlength": gridlength,
        "resol": resol,
        "step_factor": step_factor,
        "save_number": save_number,
    }
    for name, value in values.items():
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise TypeError(f"{name} must be a real number.")
        if not np.isfinite(value):
            raise ValueError(f"{name} must be finite.")

    if tmax <= 0:
        raise ValueError("tmax must be greater than zero.")
    if gridlength <= 0:
        raise ValueError("gridlength must be greater than zero.")
    if resol <= 0 or int(resol) != resol:
        raise ValueError("resol must be a positive integer.")
    if step_factor <= 0:
        raise ValueError("step_factor must be greater than zero.")
    if save_number <= 0 or int(save_number) != save_number:
        raise ValueError("save_number must be a positive integer.")

    resol = int(resol)
    save_number = int(save_number)

    # Compute spatial step size
    dx = gridlength/resol
    dt = dx**2/np.pi  # Default timestep

    # Compute the minimum number of time steps required
    min_num_steps = int(tmax/dt) + 1
    min_num_steps_adjusted = max(1, int(min_num_steps // step_factor))

    if save_number >= min_num_steps_adjusted:
        print("WARNING: The min_num_steps_int was adjusted to match save_number")
        actual_num_steps = save_number
        its_per_save = 1
    else:
        # Ensure the total number of steps is a multiple of save_number
        rem = min_num_steps_adjusted % save_number
        actual_num_steps = min_num_steps_adjusted + (save_number - rem) % save_number
        its_per_save = actual_num_steps // save_number

    # Compute final timestep
    ht = tmax / actual_num_steps
    
    return ht, its_per_save, actual_num_steps

########### Functions to check possible overlap between the objects
#############################################################################
def overlap(positionCen, rmax=5.35854, warn=0):
    """
    Checks if any solitons in the list overlap.

    Parameters:
    - positionCen: List of center positions of the soliton objects:
        Example: positionCen = [[x, y, z], [x2, y2, z2], ...]
    - rmax (float): Half of the maximum distance threshold.
    - warn (int): Counter for overlap warnings (default 0).

    Returns:
    - int: The total number of overlaps detected.
    """
    
    if any(len(sol) != 3 for sol in positionCen):  # Validate data structure
        raise ValueError("The position (center) list will be [x, y, z].")
    
    if len(positionCen) < 2:  # No comparisons can be made
        return warn

    positionCen = np.asarray(positionCen, dtype=float)  # Convert once for efficiency

    for k in range(len(positionCen) - 1):
        P_candidate = positionCen[k]
        P_solitons = positionCen[k+1:]

        overlNum = overlap_check(P_candidate, P_solitons, rmax=rmax)
        
        if overlNum:
            warn += overlNum
            print(f'WARNING: {overlNum} overlaps detected.')

    return warn


def overlap_check(P_candidate, P_solitons, rmax):
    """   
    Checks the number of solitons whose centers are within a distance of `2 * rmax` 
    from the center of a given `P_candidate`.

    Parameters:
        P_candidate (array-like [xc, yc, zc]): Center of the candidate soliton in codec units.
        P_solitons (array-like [[x0, y0, z0], [x1, y1, z1], ...]): Array of soliton center coordinates in codec units.
        rmax (float): Half of the maximum distance threshold.

    Returns:
        int: The number of solitons overlapping with `P_candidate`.
    """
    P_candidate = np.asarray(P_candidate, dtype=float)
    
    if P_candidate.shape != (3,):
        raise ValueError("P_candidate must be an array-like of length 3.")

    if len(P_solitons) == 0:
        return 0  # Handle empty input safely

    P_solitons = np.asarray(P_solitons, dtype=float)

    if P_solitons.shape[1] != 3:
        raise ValueError("P_solitons must be an array-like of shape (N, 3).")

    difference = P_candidate - P_solitons  # Broadcasting subtraction
    dist_sq = np.sum(difference**2, axis=1)  # Squared distance
    overlNum = np.sum(dist_sq <= (2 * rmax) ** 2)
    
    return int(overlNum)

########### Progress bar
#############################################################################
def progressbar(current_value, total_value, bar_length=20, progress_char='#'): 
    """
    Display a progress bar in the console.
    :param current_value: Current progress value.
    :param total_value: Total value for completion.
    :param bar_length: Length of the progress bar.
    :param progress_char: Character used to fill the progress bar.
    """
    if total_value == 0:
        print("Error: total_value cannot be 0")
        return
    
    # Calculate the percentage and progress
    percentage = int((current_value / total_value) * 100)
    progress = int((bar_length * current_value) / total_value)
    
    # Build the progress bar string
    loadbar = f"Progress: [{progress_char * progress}{'.' * (bar_length - progress)}] {percentage}%"
    
    # Print the progress bar
    print(loadbar, end='\r')
    
########### Save parameters
#############################################################################

def save_parameters(model_parameters, simulation_parameters, salva_data, comp_conserv, save_name="parameters"):
    """
    Saves simulation parameters to a text file.

    Args:
        simulation_parameters (dict): Dictionary with simulation parameters.
        salva_data (dict): Dictionary containing saving parameters (must include "address").
        comp_conserv (dict): Dictionary with complementary conserved quantities.
        name (str, optional): Name of the output file (default is "parameters").

    Returns:
        None
    """
    # Define keys to exclude
    not_save = {"profiles", "address", "plott0", "Boverlap", "info"}

    # Model data dictionary
    data_model = {}
    for name, data_list in model_parameters.items():
        data_model[name] = [
            {k: v for k, v in model.items() if k not in not_save}
            for model in data_list
        ]

    # Merge dictionaries
    dictGlobal = {**data_model, **simulation_parameters, **salva_data, **comp_conserv}

    # Get address and ensure it's a valid path
    address = dictGlobal.get("address", "./")  # Default to current directory
    address = os.path.abspath(address)  # Convert to absolute path

    # Ensure directory exists
    os.makedirs(address, exist_ok=True)

    # Construct full file path
    file_path = os.path.join(address, f"{save_name}.txt")

    # Write parameters to file
    with open(file_path, 'w', encoding="utf-8") as f:
        for param_name, value in dictGlobal.items():
            if param_name not in not_save:
                if isinstance(value, dict):  # Handle nested dictionary
                    for sub_key, sub_value in value.items():
                        f.write(f"{sub_key} >>> {sub_value}\n\n")
                else:
                    f.write(f"{param_name} >>> {value}\n\n")

    print("#" * 10 + " Parameters saved " + "#" * 10)
    print(f"File saved at: {file_path}")

    return None


## MASS VALUE
#######################
def massVal(r, sigtot, gamma=0, fac=4*np.pi, kind='quadratic', fill_value="extrapolate"):
    """
    Calculates the configuration mass based on the radial density profiles sigtot.

    Parameters:
    - r: array-like, radial positions where `sigtot` is evaluated.
    - sigtot: array-like, density profiles as a function of `r`.
    - gamma: int parameter, 1 for radial polarization, 0 for the rest (default is 0).
    - kind: string, interpolation method for `sigtot` (default is 'quadratic').

    Returns:
    - Mas: float, calculated mass after integration.

    Note:
    For multifrequency cases:
    sigtot = sigs[0]**2 + sigs[1]**2 + sigs[2]**2
    """

    # Ensure r and sigtot are sorted
    if not all(r[i] < r[i+1] for i in range(len(r)-1)):
        raise ValueError("Input array `r` must be strictly increasing.")
    
    # Ensure that gamma is either 0 or 1
    if gamma not in [0, 1]:
        raise ValueError("gamma value must be 0 or 1.")
    
    # Interpolate sigtot
    sigF = interp1d(r, sigtot, kind=kind, fill_value=fill_value)

    # Define the integrand Bf(r)
    Bf = lambda r: r**(2*(gamma+1)) * sigF(r)**2

    # Integration limits
    rmin, rmax = r[0], r[-1]

    # Perform the numerical integration
    integral_result, error = quad(Bf, rmin, rmax)
    
    if error > 1e-5:  # Large error margin, can be adjusted
        print(f"Warning: The integral may not have converged well. Error estimate: {error}")

    # Calculate the mass (Mas)
    Mas = fac*integral_result  # masa: c*hb/(G*m*Lambda^(1/2))  -> Lambda=4pi m^3/Mp^2
    return Mas


########### Recovering parameters
#############################################################################
def read_parameter(address, name):
    type_data = {
        "lambda_value": int, 
        "num_threads": int,
        "gridlength": float,
        "resol": int,
        "step_factor": float,
        "t0": float,
        "tmax": float,
        "rmax": float,
        "Plim": float,
        "cmass": float,
        "methodEnerg": int,
        "format": str,
        "save_number": int,
        "grid": bool,
        "save_rho": bool,
        "save_psi": bool,
        "save_phi": bool,
        "save_plane": bool,
        "save_energies": bool,
        "save_line": bool,
        "Numb_Part": bool,
        "Energ": bool,
        "Pi": bool,
        "Ji": bool,
        "Frequency": bool
    }

    with open(address, 'r') as file:
        for line in file:
            if '>>>' in line:
                key, value = line.strip().split('>>>')
                key = key.strip()
                value = value.strip()
                if key == name:
                    expected_type = type_data.get(name, str)
                    # Special handling for boolean values
                    if expected_type == bool:
                        return value.lower() == 'true'
                    try:
                        return expected_type(value)
                    except ValueError:
                        print(f'Error converting value for "{name}".')
                        return None
    # If not found
    print(f'Parameter "{name}" not found.')
    return None

def give_parameter(address, name):
    value = read_parameter(address, name)
    if value is not None:
        return value
    else:
        print(f'Could not retrieve the value of "{name}".')
