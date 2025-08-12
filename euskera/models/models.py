# euskera v1.0
# models file

import sys
import os

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ""))
sys.path.append(parent_dir)

import models.soliton_model as sm
import models.gaussiana_model as gm
import models.ell_model as ell_m

###################################################################################################

########### Models Class (to choosed the model)
#############################################################################
class Models():
    """ 
    Main class used to choosed the model
    """
    
    # Class attribute
    
    ###################################################################################################

    def __init__(self, info=False, **kwargs):
        """
        Initialize the Models class.

        Parameters:
        - info (bool): Whether to print initialization details.
        - kwargs (dict): Dictionary where keys are model names, and values are lists of parameter update dictionaries.
        """
        # extracting the names
        self.name = [name.lower() for name in kwargs.keys()]
        
        # Validate model/s existence and creating parameter-empty models dictionary
        self.parameters_mod = {name: dict_names(name) for name in self.name}
        
        # Creating type-data model dictionary
        models_type = {name: dict_type(name) for name in self.name}
        
        # Checking type and updating parameters
        for name in self.name:
            self.parameters_mod[name] = update_parameters(kwargs[name],  # parameters_update
                                                          self.parameters_mod[name],  # modelo_parameters
                                                          models_type[name])  # model_type
        # Instance attributes
        self.model = {
            "soliton": sm.Soli_Model(self.parameters_mod.get("soliton", {})),  # Create an instance of soliton_model
            "gaussian_function": gm.Gaussiana_Model(self.parameters_mod.get("gaussian_function", {})), # Create an instance of gaussian_model
            "ell_boson": ell_m.ell_Model(self.parameters_mod.get("ell_boson", {})) # Create an instance of gaussian_model
        }
        
        if info:
            print(f"Models initialized: {self.name} with parameters: {self.parameters_mod}")
            
            
    def call_model(self, field_components, parameters_simulation):
        """
        Calls the appropriate model and returns computed values.

        Parameters:
            field_components: number of the components of the field.
            parameters_simulation: Simulation parameters.

        Returns:
            tuple: Arrays representing computed field and its density.
        """
        grid = False
        xarray = yarray = zarray = distarray = psi = rho_i = None
    
        if "soliton" in self.name:
            model_self = self.model["soliton"]
            [xarray, yarray, zarray, distarray], [psi, rho_i] = model_self.PsiInic(field_components, parameters_simulation)
            grid = True
            
        if "gaussian_function" in self.name:
            model_self = self.model["gaussian_function"]
            if grid and psi is not None:
                [psi, rho_i] = model_self.GaussSolProf(field_components, parameters_simulation, psi, xarray, yarray, zarray)
            else:
                [xarray, yarray, zarray, distarray], [psi, rho_i] = model_self.GaussSolProf(field_components, parameters_simulation,
                                                                                            psi, xarray, yarray, zarray)
                grid = True
                
        #if "ell_boson" in self.name:
        #    model_self = self.model["ell_boson"]
        #    if grid and psi is not None:
        #        grid_data = [xarray, yarray, zarray, distarray] 
        #        [psi, rho_i] = model_self.PsiInic(field_components, parameters_simulation, grid=grid_data)
        #    else:
        #        [xarray, yarray, zarray, distarray], [psi, rho_i] = model_self.PsiInic(field_components, parameters_simulation)
        #        grid = True
            
        return [xarray, yarray, zarray, distarray], [psi, rho_i]
    
    
    
########### Extra functions (Outside of the class)
#############################################################################

def update_parameters(parameters_update, modelo_parameters, model_type):
    """
    Updates the given modelo_parameters dictionary with new values from parameters_update.
    
    Parameters:
        parameters_update (list of dicts): New values to update.
        modelo_parameters (dict): Current parameters dictionary.
        model_type (dict): Expected types for each parameter.
    
    Returns:
        dict: Updated modelo_parameters.
    """
    
    if not isinstance(parameters_update, (list, tuple)):
        raise TypeError(f"Expected a list or tuple of parameters_update, got {type(parameters_update)}.")
    
    for parameters_update_temp in parameters_update:
        if not isinstance(parameters_update_temp, dict):
            raise TypeError(f"Each item in parameters_update must be a dictionary, got {type(parameters_update_temp)}.")

        for name, element_update in parameters_update_temp.items():
            expected_type = model_type.get(name)
                
            if expected_type is None:
                raise ValueError(f"Unknown parameter: '{name}' is not a valid model parameter.")

            if not isinstance(element_update, expected_type):
                raise TypeError(f"The parameter '{name}' must be of type {expected_type}, but got {type(element_update)}.")

            modelo_parameters[name].append(element_update)  # Append new value
    
    return modelo_parameters

def dict_names(name):
    """
    Returns a dictionary of default parameters based on the given model name.
    """
    
    # Creating empty model dictionaries
    param_name_gauss = ("positions_gaussiana", "amplitude", "sigma")
    dict_gauss_func = {key: [] for key in param_name_gauss}
    
    param_name_soliton = ("profiles", "positions", "velocities", "betas", "phases", "alphas", "dr")
    dict_soliton = {key: [] for key in param_name_soliton}
    
    # Mapping model names to their respective type dictionaries
    default_parameters = {
        "soliton": dict_soliton,
        "gaussian_function": dict_gauss_func,
        "ell_boson": dict_soliton
    }
    dict_prop = default_parameters.get(name)
    
    # cheking
    if dict_prop is None:
        raise ValueError(f"Unknown model. Available models: {list(default_parameters.keys())}.")
    
    return dict_prop

def dict_type(name):
    """
    Returns a dictionary specifying the expected data types for different parameters
    based on the given model name.
    """
    
    # Using dictionary comprehension for clarity and consistency
    param_soliton = ("profiles", "positions", "velocities", "betas", "phases", "alphas", "dr")
    dict_soliton_type = {key: (tuple, list) for key in param_soliton}
    
    param_gauss = ("positions_gaussiana", "sigma")
    dict_gauss_func_type = {key: (tuple, list) for key in param_gauss}
    
    dict_gauss_func_type.update({
        "amplitude": (int, float)
    })
    
    # Mapping model names to their respective type dictionaries
    default_parameters_type = {
        "soliton": dict_soliton_type,
        "gaussian_function": dict_gauss_func_type,
        "ell_boson": dict_soliton_type
    }
    dict_prop = default_parameters_type.get(name)
    
    # checking
    if dict_prop is None:
        raise ValueError(f"Unknown model. Available models: {list(default_parameters_type.keys())}.")
    
    return dict_prop

def solitonProf(field_components, parameters_sol, simulation_parameters):
    """
    Compute the scalar field using the soliton model given a radial profile.
    """
    model_used = Models(**parameters_sol)
    (xarray, yarray, zarray, _), (_, rho_i) = model_used.call_model(field_components=field_components,
                                                                            parameters_simulation=simulation_parameters)
    
    import numexpr as ne
    rho = ne.evaluate("sum(rho_i, axis=0)")
    return  (xarray, yarray, zarray), rho