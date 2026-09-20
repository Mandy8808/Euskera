"""Canonical background-solution workflows.

The old ``background`` module remains a compatibility facade; new code should
import from these responsibility-focused modules.
"""
from .systems import *
from .profiles import *
from .asymptotics import *
from .solvers import *
from .multifrequency import *
from .fitting import *

__all__ = [
    "system", "systemMultifrequency", "systemMultFreqTot", "MatrizDXDu",
    "profilesFromSolut", "mainExt", "extend", "sParam", "Asy_uProf",
    "Asy_sProf", "Asy_sProf_v2", "extend_Met1", "extend_Met2", "shoot",
    "freq_shoot", "identify", "freq_shoot2", "identify2",
    "MultFreq_solveG2_vO", "MultFreq_solveG2", "fitting", "algebSyst",
]
