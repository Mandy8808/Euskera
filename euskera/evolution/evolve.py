"""MAIN FILE FOR EVOLUTION ROUTINES"""

import numexpr as ne
import numpy as np

# Check if pyFFTW is available
try:
    import pyfftw
    pyfftwOpt = True
except ImportError:
    pyfftw = None
    pyfftwOpt = False

from euskera.observables import simulation_conserv_quant as cq
from euskera.evolution import evolut_routines as ev
from euskera.core import grids as gd
from euskera.evolution import potential as pt
from euskera.models import models as md
from euskera.visualization import simulation_plots as pl
from euskera.io import save_data as sv
from euskera.tools import tools as to
from euskera.evolution.config import EvolutionConfig
from euskera.io.config import OutputConfig
from euskera.observables.config import DiagnosticsConfig

###################################################################################################


########### Main function
#############################################################################
def evolve(model_parameters,
           field_components=1,
           salva_data_update=None,
           simulation_parameters_update=None,
           comp_conserv_update=None,
           evolution_config=None,
           output_config=None,
           diagnostics_config=None,
           info=False):

    ######### CHECKING IF THE MODEL EXIST IN PALET OF MODEL
    #model_used = md.Models(modelo_name=model, parameters_update=model_parameters, info=info)
    model_used = md.Models(info=info, **model_parameters)

    ######### DEFAULT CONFIGURATIONS
    # Default simulation parameters
    simulation_parameters = {
                            "lambda_value": 0,
                            "num_threads": 1,
                            "gridlength": 10,
                            "resol": 128,
                            "step_factor": 1.,
                            "t0": 0,
                            "tmax": 1,
                            "rmax": 7.6,
                            "Plim": 5.6,
                            "cmass": 0,
                            "plott0": False,
                            "Boverlap": True,
                            "methodEnerg": 1,
                            "info" : False
                                  }
    # añadir check the las componentes y el numero de componentes

    if simulation_parameters_update:  # updating the default simulation parameters
        simulation_parameters = to.update_simulation_parameters(simulation_parameters_update,
                                     simulation_parameters)
    # Default save_data parameters
    salva_data = {
        "format": "npz",
        "address": "Data",
        "save_number": 10,
        "data_save" : {
            "grid": True,
            #
            "save_rho": False,
            "save_psi": False,
            "save_phi": False,
            #
            "save_plane_rho": False,
            "save_plane_psi": False,
            "save_plane_phi": False,
            #
            "save_line_rho": True,
            "save_line_psi": True,
            "save_line_phi": True,
            #
            "save_energies": True}
        }
    if salva_data_update:
        salva_data = to.update_simulation_parameters(salva_data_update, salva_data)

    comp_conserv = {
        "Numb_Part": True,
        "Energ": True,
        "Pi": False,
        "Ji": False,
        "Frequency": False
    }
    if comp_conserv_update:
        comp_conserv = to.update_simulation_parameters(comp_conserv_update, comp_conserv)

    # Legacy mappings are applied first for compatibility.  A canonical
    # dataclass is an explicit, newer choice and therefore wins on conflicts.
    if evolution_config is not None:
        if isinstance(evolution_config, dict):
            evolution_config = EvolutionConfig.from_legacy(evolution_config)
        if not isinstance(evolution_config, EvolutionConfig):
            raise TypeError("evolution_config must be an EvolutionConfig")
        simulation_parameters.update(evolution_config.to_dict())
    if output_config is not None:
        if isinstance(output_config, dict):
            output_config = OutputConfig.from_legacy(output_config)
        if not isinstance(output_config, OutputConfig):
            raise TypeError("output_config must be an OutputConfig")
        output_values = output_config.to_dict()
        output_values["data_save"] = dict(output_config.data_save)
        salva_data.update(output_values)
    if diagnostics_config is not None:
        if isinstance(diagnostics_config, dict):
            diagnostics_config = DiagnosticsConfig.from_legacy(diagnostics_config)
        if not isinstance(diagnostics_config, DiagnosticsConfig):
            raise TypeError("diagnostics_config must be a DiagnosticsConfig")
        comp_conserv.update(diagnostics_config.to_dict())

    ######################### Saving the parameters
    to.save_parameters(model_parameters, simulation_parameters, salva_data, comp_conserv, save_name="parameters")

    ######################### Check if pyFFTW is available
    if not pyfftwOpt:
        print("WARNING: pyFFTW not available, using NumPy instead.")

    ######################### Set the number of threads for NumExpr parallelization
    num_threads = simulation_parameters.get('num_threads')
    ne.set_num_threads(num_threads)

    ######################### Generating the save objects asociated to data_save
    data_save = salva_data.get("data_save")
    address = salva_data.get("address")
    formt = salva_data.get("format")
    data_save_obj = sv.data_Objgenerator(data_save=data_save, address=address, format=formt, comp_conserv=comp_conserv)

    ######################### Initialize wavefunction and density
    (xarray, yarray, zarray, distarray), (psi, rho_i) = model_used.call_model(field_components=field_components,
                                                                            parameters_simulation=simulation_parameters)
    ######################### POTENTIAL at t=0
    # Compute real/complex Fourier-space grids
    gridlength = simulation_parameters.get("gridlength")
    resol = simulation_parameters.get("resol")
    _, rkarray2 = gd.KGrid(gridlength=gridlength, resol=resol, realspace=True)
    kvec, karray2 = gd.KGrid(gridlength=gridlength, resol=resol)

    # Compute initial potential field
    cmass = simulation_parameters.get("cmass")
    (phisp, rho, obj) = pt.Upotential(field_components, rho_i, distarray, rkarray2, num_threads,
                                    cmass=cmass, resol=resol)
    # obj -> [rfft_rho, irfft_phi]

    # Optional: Plot initial density and potential
    plott0 = simulation_parameters.get("plott0")
    if plott0:
        pl.ShowPlaneProf(rho, xarray[:, 0, 0], yarray[0, :, 0], Z=None, indX=None, indY=None, indZ=None)
        pl.ShowPlaneProf(phisp, xarray[:, 0, 0], yarray[0, :, 0], Z=None, indX=None, indY=None, indZ=None)
    ################################################################################################################

    ########################## Conserved
    if comp_conserv.get("Frequency"):  # Identifying the position of the max
        max_index = [np.argmax(psi_i) for psi_i in psi]
        max_pos = [np.unravel_index(max_index_i, psi[0].shape) for max_index_i in max_index]
        if info:
            print("Index of the maximum of psi: ", max_pos)
    else:
        max_pos = None

    methodEnerg = simulation_parameters.get("methodEnerg")
    data = [psi, rho, phisp, distarray, karray2, kvec]
    cData = cq.Conserv(data, comp_conserv, simulation_parameters, obj2=None, methodEnerg=methodEnerg, max_pos=max_pos)
    if info:
        print("Initial conserved quantities:", cData)
    ################################################################################################################

    ########################## Saving
    data = ([xarray, yarray, zarray], rho, psi, phisp, cData)
    sv.fdata_save(ti=0, data=data, data_save_obj=data_save_obj, resol=resol, end=False)
    ################################################################################################################

    ######################### Compute time step parameters
    tmax = simulation_parameters.get("tmax")
    step_factor = simulation_parameters.get("step_factor")
    save_number = salva_data.get("save_number")
    ht, its_per_save, num_steps = to.dtime(tmax, gridlength, resol, step_factor, save_number)

    ################################################################################################################

    ######################### Initialize simulation fields and parameters
    lambda_value = simulation_parameters.get("lambda_value")
    halfstepornot = True  # True for a half step False for a full step
    fields = [phisp, psi, rho]
    param = [num_steps, ht, halfstepornot, its_per_save, num_threads, cmass, resol, lambda_value, max_pos]
    distDat = [distarray, karray2, rkarray2]

    # Evolve the system
    rho, phisp = ev.PKP(field_components, fields, param, distDat,
                        num_steps, obj, kvec, simulation_parameters,
                        comp_conserv, data_save_obj=data_save_obj,
                        info=info)
    ################################################################################################################
    return None
