"""MODELS REGISTRY AND HANDLER FOR EUSKERA SIMULATIONS"""

from euskera.models import ell_model as ell_m
from euskera.models import gaussiana_model as gm
from euskera.models import proca_model as proc_m
from euskera.models import soliton_model as sm

###################################################################################################

# Registry of available models
#############################################################################

MODEL_REGISTRY = {
    "soliton": sm.Soli_Model,
    "gaussian_function": gm.Gaussiana_Model,
    "ell_boson": ell_m.ell_Model,
    "proca": proc_m.proca_Model
}

########### Models Class (to choosed the model) #############################
#############################################################################
class Models():
    """ 
    Main class used to composable physical models

    All models must implement:

        apply(field_components, parameters_simulation,
              grid_data=None, psi=None)

    and return:

        grid_data, psi, rho_i
    """
    
    # Class attribute
    ###################

    def __init__(self, info=False, **kwargs):
        """
        Initialize the Models class.

        Parameters:
        - info (bool): Whether to print initialization details.
        - kwargs (dict): Dictionary where keys are model names, and values are lists of parameter update dictionaries.
        """
        if not kwargs: raise ValueError("At least one model must be specified.")

        # extracting the names
        self.name = [name.lower() for name in kwargs.keys()]
        
        # Creating type-data model dictionary
        models_type = {name: dict_type(name) for name in self.name}
        
        # Checking type and updating parameters
        self.parameters_mod = {} 
        for name in self.name:
            self.parameters_mod[name] = update_parameters(kwargs[name],  # parameters_update
                                                          models_type[name])  # model_type
        # Instance attributes
        self.model = {
            name: MODEL_REGISTRY[name](self.parameters_mod.get(name, {}))
            for name in self.name
        }

        if info:
            print(f"Models initialized in order: {self.name}")
            print(f"Parameters: {self.parameters_mod}")
            
    def call_model(self, field_components, parameters_simulation):
        """
        Calls the appropriate model sequentially applies models in user-defined order..            

        Parameters:
            field_components: number of the components of the field.
            parameters_simulation: Simulation parameters.

        Returns:
            grid_data, [psi, rho_i]
        """

        grid_data = None
        psi = None
        rho_i = None

        for name in self.name:
            model_self = self.model[name]

            grid_data, psi, rho_i = model_self.apply(
                field_components,
                parameters_simulation,
                grid_data=grid_data,
                psi=psi
            )
        
        return grid_data, [psi, rho_i]
    
    
########### Extra functions (Outside of the class)
#############################################################################

def update_parameters(parameters_update, model_type):
    """
    Give the modelo_parameters dictionary from parameters_update.
    
    Parameters:
        parameters_update (list of dicts): New values to update.
        model_type (dict): Expected types for each parameter.
    
    Returns:
        dict: Updated modelo_parameters.
    """
    
    if not isinstance(parameters_update, (list, tuple)): raise TypeError(f"Expected a list or tuple of parameters_update, got {type(parameters_update)}.")
    
    modelo_parameters = {}
    for parameters_update_temp in parameters_update:
        if not isinstance(parameters_update_temp, dict):
            raise TypeError(f"Each item in parameters_update must be a dictionary, got {type(parameters_update_temp)}.")

        for name, element_update in parameters_update_temp.items():
            expected_type = model_type.get(name)
                
            if expected_type is None: raise ValueError(f"Unknown parameter: '{name}' is not a valid model parameter.")
            if not isinstance(element_update, expected_type): raise TypeError(f"The parameter '{name}' must be of type {expected_type}, but got {type(element_update)}.")

            modelo_parameters.setdefault(name, []).append(element_update)
    
    return modelo_parameters

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

    dict_ell_boson_type = {key: (tuple, list) for key in param_soliton}
    dict_ell_boson_type.update({
        "ell": (tuple, list)
    })

    dict_proca_type = {key: (tuple, list) for key in param_soliton}
    dict_proca_type.update({
        "polarization": (tuple, list)
    })
    
    # Mapping model names to their respective type dictionaries
    default_parameters_type = {
        "soliton": dict_soliton_type,
        "gaussian_function": dict_gauss_func_type,
        "ell_boson": dict_ell_boson_type,
        "proca": dict_proca_type
    }
    dict_prop = default_parameters_type.get(name)
    
    # checking
    if dict_prop is None: raise ValueError(f"Unknown model. Available models: {list(default_parameters_type.keys())}.")
    return dict_prop

def solitonProf(field_components, parameters_sol, simulation_parameters):
    """
    Compute the scalar field using the soliton model given a radial profile.
    """
    model_used = Models(**parameters_sol)
    grid_data, (_, rho_i) = model_used.call_model(field_components=field_components,
                                                                    parameters_simulation=simulation_parameters)
    
    xarray, yarray, zarray, _ = grid_data
    
    import numexpr as ne
    rho = ne.evaluate("sum(rho_i, axis=0)")
    
    return  (xarray, yarray, zarray), rho