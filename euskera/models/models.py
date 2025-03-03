# euskera v1.0
# models file

import sys
import os

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ""))
sys.path.append(parent_dir)

import models.soliton_model as sm
###################################################################################################

########### Models Class (to choosed the model)
#############################################################################
class Models():
    """ 
    Main class used to choosed the model
    """
    
    # Class attribute
    dict_soliton_type = {
        "profiles": (tuple, list), "positions": (tuple, list), "velocities": (tuple, list),
        "betas":(tuple, list), "phases": (tuple, list), "alphas": (tuple, list), "dr": (tuple, list)
    }
    dict_gauss_func_type = {}  # complete with the gauss parameters
    
    default_parameters_type = {
        "soliton": dict_soliton_type,
        "gaussian_function": dict_gauss_func_type}
    
    ###################################################################################################
    
    def __init__(self, modelo_name, parameters_update, info=False):
        """
        Initialize the storage class.

        Parameters:
        - modelo_name (str): Base name used to identify the model.
        - parameters_update (list of dict): List of parameter update dictionaries.
        - info (bool): Print information about file saving.
        """
        
        # Creating empty model dictionary
        self.dict_gauss_func = {}
        self.dict_soliton = {"profiles": [], "positions": [], "velocities": [], 
                    "betas": [], "phases": [], "alphas": [], "dr": []}
        
        self.default_parameters = {"soliton": self.dict_soliton,
                          "gaussian_function": self.dict_gauss_func}
        
        # Starting the updating 
        modelo_name = modelo_name.lower()
        # Validate model existence
        modelo_parameters = self.default_parameters.get(modelo_name)
        if modelo_parameters is None:
            raise ValueError(f"Unknown model. Available models: {list(self.default_parameters.keys())}.")
        modelo_type = self.default_parameters_type[modelo_name.lower()]
        
        # Update parameters
        for parameters_update_temp in parameters_update:
            for name, element_update in parameters_update_temp.items():
                expected_type = modelo_type.get(name)
                
                if expected_type is None:
                    raise ValueError(f"Unknown parameter: '{name}' is not a valid model parameter.")

                if not isinstance(element_update, expected_type):
                    raise ValueError(f"The parameter '{name}' must be of type {expected_type}.")
                
                modelo_parameters[name].append(element_update)  # add to the empty list
        
        # Instance attributes
        self.name = modelo_name.lower()
        self.parameters_mod = modelo_parameters
        self.soli_model = sm.Soli_Modelo(self.parameters_mod)  # Create an instance of Soli_Modelo
        
        if info:
            print(f"The model used is: {modelo_name.lower()}, with parameters: {modelo_parameters}")
            
        
    def call_model(self, field_components, simulation_parameters):
        """ 
        """
        if self.name == "soliton":
            
            [xarray, yarray, zarray, distarray], [psi, rho_i] = self.soli_model.PsiInic(field_components, simulation_parameters)
        else:
            pass
            
        return [xarray, yarray, zarray, distarray], [psi, rho_i]
            