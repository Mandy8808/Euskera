"""Public API for Euskera's numerical evolution workflows."""

__version__ = "1.1.0"

# Main evolution routines
from .evolution.evolve import evolve
from .evolution.evolut_routines import PKP
from .core.grids import RealGrid, KGrid, meshgrid
from .evolution.potential import Upotential
from .observables.simulation_conserv_quant import Conserv, Npar, Energ, centpotetE, selfinterCondensateE, selfinterFieldE, kintE, Pi
from .observables.frequency import main_frequency, frequMet1, frequMet2

# Video visualization (optional import)
from .visualization.video_make import Visualization

# Tools
from .tools.tools import update_simulation_parameters, dtime, overlap, overlap_check, progressbar, massVal, read_parameter, give_parameter

# Data saving
from .io.save_data import data_Objgenerator, fdata_save, nameData, JoinFilesInOneZip, StoreSolution

# Plot configurations
from .visualization.configuration import general, FigParam, LineParam, axesParam, labelParam, legendParam, fontParam
from .visualization.plot_tools import colorBar_and_normaliz, colored_line
from .visualization.simulation_plots import ShowPlaneProf, PlaneProf, Plot3DCorrProf, imagshow2D
from .visualization.background_plots import plotUsingPerf, plotUsingDiscSol, plotUsingSol, plotPerf

# Models
from .models.models import Models, update_parameters, dict_type, solitonProf
from .models.soliton_model import Soli_Model, build_soliton, initsoliton_kernel
from .models.ell_model import ell_Model, build_ell, compute_sph_harm_grid, initell_kernel
from .models.gaussiana_model import Gaussiana_Model, build_1d_gaussian, add_gaussian
from .models.proca_model import proca_Model, build_proca, compute_polarization_vec, initproca_kernel

# Workflow-oriented canonical namespaces.
from . import backgrounds, core, evolution, io, observables, spectral, visualization
from .evolution.config import EvolutionConfig
from .io.config import OutputConfig
from .observables.config import DiagnosticsConfig

# Define available imports
__all__ = [
    "__version__",
    # Core evolution routines
    'evolve', 'PKP',
    'RealGrid', 'KGrid', 'meshgrid',
    'Upotential',
    'Conserv', 'Npar', 'Energ', 'centpotetE', 'selfinterCondensateE', 'selfinterFieldE', 'kintE', 'Pi',
    'main_frequency', 'frequMet1', 'frequMet2',
    
    # Video visualization
    'Visualization',
    
    # Tools
    'update_simulation_parameters', 'dtime', 'overlap', 'overlap_check', 'progressbar', 'massVal', 'read_parameter',
    'give_parameter',
    
    # Data saving
    'data_Objgenerator', 'fdata_save', 'nameData', 'JoinFilesInOneZip', 'StoreSolution',
    
    # Plotting
    'general', 'FigParam', 'LineParam', 'axesParam', 'labelParam', 'legendParam', 'fontParam',
    'colorBar_and_normaliz', 'ShowPlaneProf', 'PlaneProf', 'Plot3DCorrProf', 'imagshow2D', 'colored_line',
    'plotUsingPerf', 'plotUsingDiscSol', 'plotUsingSol', 'plotPerf',
    
    # main
    'Models', 'update_parameters', 'dict_type', 'solitonProf',
    
    # soliton
    'Soli_Model', 'build_soliton', 'initsoliton_kernel',

    # ell_boson
    'ell_Model', 'build_ell', 'compute_sph_harm_grid', 'initell_kernel', 
    
    # guassiana
    'Gaussiana_Model', 'build_1d_gaussian', 'add_gaussian',
    
    # proca
    'proca_Model', 'build_proca', 'compute_polarization_vec', 'initproca_kernel',

    # Workflow namespaces
    'backgrounds', 'core', 'evolution', 'io', 'observables',
    'visualization',
    'EvolutionConfig', 'OutputConfig', 'DiagnosticsConfig',
]
from euskera.io.schedule import All, Last, TimeRange, Final, SaveRule
from euskera.io.selected_output import read_output
__all__ += ["All", "Last", "TimeRange", "Final", "SaveRule", "read_output"]
