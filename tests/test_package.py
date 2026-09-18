import importlib

import numpy as np
import pytest

import euskera
from euskera.models.models import Models


def test_package_import_exposes_public_api():
    assert callable(euskera.evolve)
    assert callable(euskera.dtime)
    assert Models is euskera.Models


def test_gaussian_model_initializes_small_wavefunction():
    model = Models(
        gaussian_function=[
            {
                "positions_gaussiana": [0.0, 0.0, 0.0],
                "amplitude": 1.0,
                "sigma": [1.0, 1.0, 1.0],
            }
        ]
    )

    grid_data, (psi, rho_i) = model.call_model(
        field_components=1,
        parameters_simulation={
            "gridlength": 4.0,
            "resol": 4,
            "lambda_value": 0.0,
            "num_threads": 1,
            "Boverlap": True,
            "Plim": 5.6,
            "rmax": 7.6,
        },
    )

    assert len(grid_data) == 4
    assert psi.shape == (1, 4, 4, 4)
    assert rho_i.shape == (1, 4, 4, 4)
    assert np.isfinite(psi).all()
    assert np.isfinite(rho_i).all()


@pytest.mark.parametrize("module_name", ["euskera.main.main", "euskera.models.models"])
def test_internal_modules_import_as_package(module_name):
    assert importlib.import_module(module_name)
