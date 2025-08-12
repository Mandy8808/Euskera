# modulo/main/__init__.py

# Import specific functions and classes from submodules
from .main import evolve
from .evolut_routines import PKP
from .grids import RealGrid, KGrid, meshgrid
from .potential import Upotential
from .conserv_quant import Conserv, Npar, Energ, centpotetE, selfinterCondensateE, selfinterFieldE, kintE, Pi
from .frequency import main_frequency, frequMet1, frequMet2

# Explicitly define what should be imported when using `from modulo.main import *`
__all__ = [
           # main 
           'evolve', 'PKP', 
           
           # grids
           'RealGrid', 'KGrid', 'meshgrid',
           
           # initial configuration
           'Upotential',
           
           # conserved
           'Conserv', 'Npar', 'Energ, ''centpotetE', 'selfinterCondensateE', 'selfinterFieldE', 'kintE', 'Pi',
           
           # frequency
           'main_frequency', 'frequMet1', 'frequMet2'
           
           ]

