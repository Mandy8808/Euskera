"""Checks for the standalone spectral package layout."""

import importlib.util

from euskera.spectral import backgroundOper, cheb, spectrum
from euskera.visualization.spectral_plot import plotImag


def test_spectral_namespace_exports_core_functions():
    assert callable(cheb)
    assert callable(backgroundOper)
    assert callable(spectrum)
    assert callable(plotImag)
    assert importlib.util.find_spec("euskera.numerics") is None
