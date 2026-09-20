"""Legacy compatibility namespace; use core, evolution and observables."""
from euskera.evolution import evolve, PKP, Upotential
from euskera.core import RealGrid, KGrid, meshgrid
from euskera.observables import Conserv, Npar, Energ, centpotetE, selfinterCondensateE, selfinterFieldE, kintE, Pi, main_frequency, frequMet1, frequMet2
__all__ = ["evolve", "PKP", "RealGrid", "KGrid", "meshgrid", "Upotential", "Conserv", "Npar", "Energ", "centpotetE", "selfinterCondensateE", "selfinterFieldE", "kintE", "Pi", "main_frequency", "frequMet1", "frequMet2"]
# Module aliases keep monkeypatching and introspection compatible with legacy imports.
from euskera.evolution import potential, evolut_routines
from euskera.core import grids
from euskera.observables import conserv_quant, frequency
