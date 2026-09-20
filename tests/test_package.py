import importlib

import importlib.util

import numpy as np
import pytest

import euskera
import euskera.backgrounds as backgrounds
import euskera.evolution as evolution
from euskera.evolution import potential
from euskera.observables.simulation_conserv_quant import Npar
from euskera.core.grids import KGrid, RealGrid
from euskera.models.models import Models
from euskera.evolution import EvolutionConfig
from euskera.io import OutputConfig
from euskera.observables import DiagnosticsConfig
from euskera.backgrounds.profiles import profilesFromSolut


def test_package_import_exposes_public_api():
    assert callable(euskera.evolve)
    assert callable(euskera.dtime)
    assert Models is euskera.Models


def test_background_package_imports():
    assert backgrounds.__name__ == "euskera.backgrounds"
    assert callable(backgrounds.system)


@pytest.mark.parametrize("module", [euskera, evolution, backgrounds])
def test_public_exports_are_defined(module):
    missing = [name for name in module.__all__ if not hasattr(module, name)]
    assert not missing, f"{module.__name__} has missing exports: {missing}"

    namespace = {}
    exec(f"from {module.__name__} import *", namespace)
    assert all(name in namespace for name in module.__all__)


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


def test_potential_is_finite_and_has_zero_mean_for_localized_density(monkeypatch):
    monkeypatch.setattr(potential, "pyfftw", None)
    monkeypatch.setattr(potential, "pyfftwOpt", False)

    (_, _, _,), distarray = RealGrid(gridlength=4.0, resol=8)
    _, rkarray2 = KGrid(gridlength=4.0, resol=8, realspace=True)
    coordinates = np.linspace(-1.75, 1.75, 8)
    x, y, z = np.meshgrid(coordinates, coordinates, coordinates, indexing="ij")
    rho_i = np.exp(-(x**2 + y**2 + z**2))[None, ...]

    phisp, rho, _ = potential.Upotential(
        field_components=1,
        rho_i=rho_i,
        distarray=distarray,
        rkarray2=rkarray2,
        num_threads=1,
        resol=8,
    )

    assert np.isfinite(phisp).all()
    assert np.isfinite(rho).all()
    assert np.isclose(phisp.mean(), 0.0, atol=1e-12)


def test_particle_number_matches_discrete_norm():
    psi = np.array([[[[1.0 + 0.0j, 2.0 + 0.0j]]]])
    total, components = Npar(psi, Vcell=0.5)

    assert components == [2.5]
    assert total == 2.5


@pytest.mark.parametrize(
    ("arguments", "exception"),
    [
        ((0.0, 10.0, 128, 1.0, 10), ValueError),
        ((1.0, 0.0, 128, 1.0, 10), ValueError),
        ((1.0, 10.0, 0, 1.0, 10), ValueError),
        ((1.0, 10.0, 128, 0.0, 10), ValueError),
        ((1.0, 10.0, 128, 1.0, 0), ValueError),
        ((1.0, 10.0, 128.5, 1.0, 10), ValueError),
    ],
)
def test_dtime_rejects_invalid_arguments(arguments, exception):
    with pytest.raises(exception):
        euskera.dtime(*arguments)


def test_dtime_returns_integer_save_interval():
    ht, its_per_save, num_steps = euskera.dtime(1.0, 10.0, 16, 1.0, 4)

    assert ht > 0
    assert its_per_save == int(its_per_save)
    assert num_steps % 4 == 0


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


def test_typed_configurations_validate_and_round_trip_legacy():
    evolution = EvolutionConfig(gridlength=4.0, resol=8, tmax=0.1)
    output = OutputConfig(address="runs/test", save_number=2)
    diagnostics = DiagnosticsConfig(Pi=True)

    assert EvolutionConfig.from_legacy(evolution.to_dict()) == evolution
    assert OutputConfig.from_legacy(output.to_dict()) == output
    assert DiagnosticsConfig.from_legacy(diagnostics.to_dict()) == diagnostics
    with pytest.raises(ValueError):
        EvolutionConfig.from_legacy({"not_a_parameter": 1})


def test_typed_configurations_are_public_canonical_imports():
    assert euskera.EvolutionConfig is EvolutionConfig
    assert euskera.OutputConfig is OutputConfig
    assert euskera.DiagnosticsConfig is DiagnosticsConfig


def test_observables_do_not_expose_redundant_facade_modules():
    assert importlib.util.find_spec("euskera.observables.mass") is None
    assert importlib.util.find_spec("euskera.observables.energy") is None


def test_evolution_does_not_expose_redundant_split_step_module():
    assert importlib.util.find_spec("euskera.evolution.split_step") is None


def test_background_profile_is_deterministic_and_finite():
    result = profilesFromSolut(
        (0.5, 1.0, 0, 0.0, [], None, "DOP853", 1e-8, 1e-9, 0.0),
        Nptos=8,
    )
    energy, mass, radius, profiles_f = result[:4]
    assert radius.shape == (8,)
    assert np.isfinite(radius).all()
    assert np.isfinite(profiles_f[0]).all()
    assert np.isfinite(energy).all()
    assert np.isfinite(mass).all()


@pytest.mark.parametrize("module_name", ["euskera.evolution.evolve", "euskera.models.models"])
def test_internal_modules_import_as_package(module_name):
    assert importlib.import_module(module_name)
