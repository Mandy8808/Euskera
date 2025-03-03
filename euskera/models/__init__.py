# modulo/models/__init__.py

from .models import Models
from .soliton_model import Soli_Modelo, InitialSolitonsProf, initsolitonInt, initsoliton
#from .modelo_gaussiana import gaussian, initGaussianInt, GaussSolProf


__all__ = [
    # main
    'Models',
    
    # soliton
    'Soli_Modelo', 'InitialSolitonsProf', 'initsolitonInt', 'initsoliton'
    
    ]