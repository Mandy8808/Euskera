"""Small deterministic checks for the public scientific workflows."""

import numpy as np

import euskera
from euskera.backgrounds import profilesFromSolut, system, systemMultifrequency
from euskera.backgrounds.systems import system as systems_system
from euskera.backgrounds.shooting import freq_shoot, shoot
from euskera.evolution.config import EvolutionConfig
from euskera.io.config import OutputConfig
from euskera.models.models import Models
from euskera.observables.conserv_quant import Npar
from euskera.observables.energy_mass import energEng, massVal


def _gaussian_model():
    return {
        "gaussian_function": [
            {
                "positions_gaussiana": [0.0, 0.0, 0.0],
                "amplitude": 1.0,
                "sigma": [1.0, 1.0, 1.0],
            }
        ]
    }


def test_canonical_evolution_finishes_and_writes_to_tmp_path(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    euskera.evolve(
        _gaussian_model(),
        evolution_config=EvolutionConfig(gridlength=4.0, resol=4, tmax=0.01),
        output_config=OutputConfig(
            address=str(tmp_path / "run"),
            save_number=1,
            data_save={
                "grid": True,
                "save_line_rho": True,
                "save_line_psi": False,
                "save_line_phi": False,
                "save_energies": False,
            },
        ),
        diagnostics_config={"Numb_Part": False, "Energ": False},
    )
    assert (tmp_path / "run" / "parameters.txt").is_file()
    outputs = list((tmp_path / "run").glob("end_*.npz"))
    assert outputs
    assert all(np.load(path).files for path in outputs)


def test_potential_density_and_particle_number_are_finite(monkeypatch):
    from euskera.evolution import potential
    from euskera.core.grids import KGrid, RealGrid

    monkeypatch.setattr(potential, "pyfftw", None)
    monkeypatch.setattr(potential, "pyfftwOpt", False)
    (_, _, _,), distances = RealGrid(gridlength=4.0, resol=4)
    _, reciprocal = KGrid(gridlength=4.0, resol=4, realspace=True)
    rho_i = np.ones((1, 4, 4, 4), dtype=float)
    phi, rho, handles = potential.Upotential(
        1, rho_i, distances, reciprocal, num_threads=1, resol=4
    )
    assert handles == (None, None)
    assert np.isfinite(phi).all() and np.isfinite(rho).all()
    total, components = Npar(rho_i.astype(complex), Vcell=0.5)
    assert total == components[0] == 32.0


def test_background_rhs_and_short_profile_observables_are_finite():
    y = np.array([0.5, 0.1, 0.2, -0.05])
    rhs = np.asarray(system(0.3, y, [0.7, 0.2, 0]))
    multi_rhs = systemMultifrequency(0.3, np.tile(y, 2), [2, 0.2, 0])
    assert rhs.shape == (4,) and np.isfinite(rhs).all()
    assert multi_rhs.shape == (8,) and np.isfinite(multi_rhs).all()
    assert shoot(0.0, 2.0) == 1.0
    interval, event_radius = freq_shoot(
        (np.array([0.2, 0.8]), np.array([0.4])), 1, 0.5, (0.0, 1.0), 0.8
    )
    assert interval == [0.0, 0.5] and event_radius == 0.8

    result = profilesFromSolut(
        (0.5, 1.0, 0, 0.0, [], None, "RK45", 1e-6, 1e-7, 0.0),
        Nptos=8,
    )
    energy, mass, radius, profiles_f = result[:4]
    assert radius.shape == (8,)
    assert np.isfinite(radius).all() and np.isfinite(profiles_f[0]).all()
    assert np.isfinite(energy).all() and np.isfinite(mass).all()


def test_canonical_background_functions_are_stable():
    y = np.array([0.5, 0.1, 0.2, -0.05])
    assert system is systems_system
    assert np.isfinite(system(0.3, y, [0.7, 0.2, 0])).all()
    r = np.linspace(0.0, 1.0, 8)
    density = np.exp(-r**2)
    assert np.isfinite(energEng(r, density, 0.2))
    assert np.isfinite(massVal(r, density))
    assert Models is euskera.Models
