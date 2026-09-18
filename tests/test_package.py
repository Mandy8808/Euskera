import importlib

import numpy as np
import pytest

import euskera
import background
from euskera.main import potential
from euskera.models.models import Models


def test_package_import_exposes_public_api():
    assert callable(euskera.evolve)
    assert callable(euskera.dtime)
    assert Models is euskera.Models


def test_background_package_imports():
    assert background.__name__ == "background"
    assert callable(background.system)


def test_potential_works_without_pyfftw(monkeypatch):
    monkeypatch.setattr(potential, "pyfftw", None)
    monkeypatch.setattr(potential, "pyfftwOpt", False)

    rho_i = np.ones((1, 4, 4, 4))
    distarray = np.ones((4, 4, 4))
    rkarray2 = np.ones((4, 4, 3))

    phisp, rho, fft_objects = potential.Upotential(
        field_components=1,
        rho_i=rho_i,
        distarray=distarray,
        rkarray2=rkarray2,
        num_threads=1,
        resol=4,
    )

    assert phisp.shape == (4, 4, 4)
    assert rho.shape == (4, 4, 4)
    assert fft_objects == (None, None)


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
