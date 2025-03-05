# modulo/models/__init__.py

from .models import Models
from .soliton_model import Soli_Model, InitialSolitonsProf, initsolitonInt, initsoliton
from .gaussiana_model import Gaussiana_Model, gaussian, initGaussianInt
#from .modelo_gaussiana import gaussian, initGaussianInt, GaussSolProf


__all__ = [
    # main
    'Models',
    
    # soliton
    'Soli_Modelo', 'InitialSolitonsProf', 'initsolitonInt', 'initsoliton',
    
    # guassiana
    'Gaussiana_Model', 'gaussian', 'initGaussianInt'
    ]