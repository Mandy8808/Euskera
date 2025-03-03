# modulo/__init__.py

# Main evolution routines
from .main.main import evolve
from .main.evolut_routines import PKP
from .main.grids import RealGrid, KGrid, meshgrid
from .main.potential import Upotential

# Video visualization (optional import)
from .video.video_make import Visualization

# Tools
from .tools.tools import update_simulation_parameters, dtime, overlap, overlap_check, progressbar

# Data saving
from .save.save_data import data_Objgenerator, fdata_save, nameData, JoinFilesInOneZip, StoreSolution

# Plot configurations
from .plots.configuration import general, FigParam, LineParam, axesParam, labelParam, legendParam, fontParam
from .plots.plot_tools import colorBar_and_normaliz, ShowPlaneProf, PlaneProf, Plot3DCorrProf, imagshow2D, colored_line

# Models
from .models.models import Models
from .models.soliton_model import Soli_Modelo, InitialSolitonsProf, initsolitonInt, initsoliton
#from .models.modelo_gaussiana import gaussian, initGaussianInt, GaussSolProf

# Define available imports
__all__ = [
    # Core evolution routines
    'evolve', 'PKP',
    'RealGrid', 'KGrid', 'meshgrid',
    'Upotential',
    
    # Video visualization
    'Visualization',
    
    # Tools
    'update_simulation_parameters', 'dtime', 'overlap', 'overlap_check', 'progressbar',
    
    # Data saving
    'data_Objgenerator', 'fdata_save', 'nameData', 'JoinFilesInOneZip', 'StoreSolution',
    
    # Plotting
    'general', 'FigParam', 'LineParam', 'axesParam', 'labelParam', 'legendParam', 'fontParam',
    'colorBar_and_normaliz', 'ShowPlaneProf', 'PlaneProf', 'Plot3DCorrProf', 'imagshow2D', 'colored_line'
    
    # Models
    'Models', 'Soli_Modelo', 'InitialSolitonsProf', 'initsolitonInt', 'initsoliton',
    # 'gaussian', 'initGaussianInt', 'GaussSolProf'
]