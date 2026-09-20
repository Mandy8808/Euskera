"""Analytic plane-wave checks through the conserved-quantity pipeline."""

import numpy as np
import pytest

from euskera.core.grids import KGrid, RealGrid
from euskera.observables.simulation_conserv_quant import Conserv


@pytest.mark.parametrize("wavevector", [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, -2, 3)])
@pytest.mark.parametrize("amplitudes", [(1.,), (1., 0.5)])
def test_kinetic_methods_match_plane_wave_energy(wavevector, amplitudes):
    length, resolution = 2 * np.pi, 8
    coordinates, distances = RealGrid(length, resolution)
    kvec, k_squared = KGrid(length, resolution)
    phase = sum(k * x for k, x in zip(wavevector, coordinates))
    wave = np.broadcast_to(np.exp(1j * phase), (resolution,) * 3)
    psi = np.stack([amplitude * wave for amplitude in amplitudes])
    rho = np.sum(np.abs(psi) ** 2, axis=0)
    # Zero potential and couplings isolate kinetic energy in Conserv/Energ.
    data = [psi, rho, np.zeros_like(rho), distances, k_squared, kvec]
    parameters = {
        "resol": resolution, "gridlength": length, "num_threads": 1,
        "cmass": 0., "lambda_value": 0.,
    }
    expected = 0.5 * sum(k * k for k in wavevector) * length**3 * sum(
        amplitude**2 for amplitude in amplitudes
    )
    energies = [
        float(Conserv(data, {"Energ": True}, parameters, methodEnerg=method)[0])
        for method in (1, 2)
    ]
    np.testing.assert_allclose(energies, expected, rtol=1e-12, atol=1e-12)
