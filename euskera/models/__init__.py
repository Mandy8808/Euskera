# modulo/models/__init__.py

from .models import Models, solitonProf
from .soliton_model import Soli_Model, InitialSolitonsProf, initsolitonInt, initsoliton
from .ell_model import ell_Model
from .gaussiana_model import Gaussiana_Model, gaussian, initGaussianInt
#from .modelo_gaussiana import gaussian, initGaussianInt, GaussSolProf


__all__ = [
    # main
    'Models', 'solitonProf',
    
    # soliton
    'Soli_Model', 'InitialSolitonsProf', 'initsolitonInt', 'initsoliton',

    # ell_boson
    'ell_Model',
    
    # guassiana
    'Gaussiana_Model', 'gaussian', 'initGaussianInt'
    ]