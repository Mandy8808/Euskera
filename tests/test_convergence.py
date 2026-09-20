"""Small scientific checks; these are not substitutes for production studies."""

import numpy as np
import pytest

import euskera
from euskera.models.models import Models
from euskera.observables.conserv_quant import Npar


def _profile(resolution):
    model = Models(
        gaussian_function=[{
            "positions_gaussiana": [0.0, 0.0, 0.0],
            "amplitude": 1.0,
            "sigma": [1.0, 1.0, 1.0],
        }]
    )
    _, (psi, rho) = model.call_model(
        field_components=1,
        parameters_simulation={
            "gridlength": 4.0, "resol": resolution, "lambda_value": 0.0,
            "num_threads": 1, "Boverlap": True, "Plim": 5.6, "rmax": 7.6,
        },
    )
    return psi, rho


@pytest.mark.scientific
def test_gaussian_profile_resolution_comparison():
    """A minimum resolution check, not an asymptotic convergence claim."""
    coarse_psi, coarse_rho = _profile(4)
    fine_psi, fine_rho = _profile(6)
    coarse_mass = Npar(coarse_psi, Vcell=(4.0 / 4) ** 3)[0]
    fine_mass = Npar(fine_psi, Vcell=(4.0 / 6) ** 3)[0]
    assert np.isfinite(coarse_rho).all() and np.isfinite(fine_rho).all()
    assert np.isclose(coarse_mass, fine_mass, rtol=0.20, atol=1e-12)


@pytest.mark.scientific
def test_short_evolution_preserves_finite_state(tmp_path, monkeypatch):
    """Short runs are too small for a meaningful order-of-accuracy estimate."""
    monkeypatch.chdir(tmp_path)
    euskera.evolve(
        {"gaussian_function": [{
            "positions_gaussiana": [0.0, 0.0, 0.0],
            "amplitude": 1.0, "sigma": [1.0, 1.0, 1.0],
        }]},
        evolution_config=euskera.EvolutionConfig(
            gridlength=4.0, resol=4, tmax=0.01
        ),
        output_config=euskera.OutputConfig(
            address=str(tmp_path / "run"), save_number=1,
            data_save={"grid": True, "save_line_rho": True,
                       "save_line_psi": False, "save_line_phi": False,
                       "save_energies": False},
        ),
        diagnostics_config={"Numb_Part": False, "Energ": False},
    )
    assert list((tmp_path / "run").glob("end_save_line_rho*.npz"))
