# euskera v1.0
# tools file

import numpy as np
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
        "format": str,
        "address": str,
        "data_save": dict,
        "Numb_Part": bool,
        "Energ": bool,
        "Pi": bool,
        "Ji": bool,
        "methodEnerg": int
    }

    for name, element_update in simulation_parameters_update.items():
        expected_type = simulation_type.get(name)

        if expected_type is None:
            raise ValueError(f"Unknown parameter: '{name}' is not a valid simulation parameter.")

        if not isinstance(element_update, expected_type):
            raise ValueError(f"The parameter '{name}' must be of type {expected_type}.")

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
    
    # Compute spatial step size
    dx = gridlength/resol
    dt = dx**2/np.pi  # Default timestep

    # Compute the minimum number of time steps required
    min_num_steps = int(tmax/dt) + 1
    min_num_steps_adjusted = min_num_steps // step_factor  # Apply step factor adjustment

    if save_number >= min_num_steps_adjusted:
        print("WARNING: The min_num_steps_int was adjusted to match save_number")
        actual_num_steps = save_number
        its_per_save = 1
    else:
        # Ensure the total number of steps is a multiple of save_number
        rem = min_num_steps_adjusted % save_number
        actual_num_steps = min_num_steps_adjusted + save_number - rem
        its_per_save = actual_num_steps / save_number  #  //

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