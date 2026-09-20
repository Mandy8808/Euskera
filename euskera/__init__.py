# modulo/__init__.py

# Main evolution routines
from .main.main import evolve
from .main.evolut_routines import PKP
from .main.grids import RealGrid, KGrid, meshgrid
from .main.potential import Upotential
from .main.conserv_quant import Conserv, Npar, Energ, centpotetE, selfinterCondensateE, selfinterFieldE, kintE, Pi
from .main.frequency import main_frequency, frequMet1, frequMet2

# Video visualization (optional import)
from .video.video_make import Visualization

# Tools
from .tools.tools import update_simulation_parameters, dtime, overlap, overlap_check, progressbar, massVal, read_parameter, give_parameter

# Data saving
from .save.save_data import data_Objgenerator, fdata_save, nameData, JoinFilesInOneZip, StoreSolution

# Plot configurations
from .plots.configuration import general, FigParam, LineParam, axesParam, labelParam, legendParam, fontParam
from .plots.plot_tools import colorBar_and_normaliz, ShowPlaneProf, PlaneProf, Plot3DCorrProf, imagshow2D, colored_line

# Models
from .models.models import Models, update_parameters, dict_type, solitonProf
from .models.soliton_model import Soli_Model, build_soliton, initsoliton_kernel
from .models.ell_model import ell_Model, build_ell, compute_sph_harm_grid, initell_kernel
from .models.gaussiana_model import Gaussiana_Model, build_1d_gaussian, add_gaussian
from .models.proca_model import proca_Model, build_proca, compute_polarization_vec, initproca_kernel

# Workflow-oriented namespaces. The legacy modules remain available for
# backwards compatibility while new code can use these stable boundaries.
from . import backgrounds, core, evolution, io, numerics, observables, visualization

# Define available imports
__all__ = [
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
    'backgrounds', 'core', 'evolution', 'io', 'numerics', 'observables',
    'visualization',
]