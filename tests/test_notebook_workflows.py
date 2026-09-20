"""Fast, deterministic equivalents of the example notebook workflows.

The example notebooks intentionally remain data-producing documents and are
not executed in CI.  These tests exercise their small, public API building
blocks instead.
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import euskera
from euskera.io import OutputConfig
from euskera.models.models import Models
from euskera.visualization.plot_tools import imagshow2D


def _gaussian_parameters():
    return {
        "gaussian_function": [
            {
                "positions_gaussiana": [0.0, 0.0, 0.0],
                "amplitude": 1.0,
                "sigma": [1.0, 1.0, 1.0],
            }
        ]
    }


def test_initial_profile_smoke_is_finite():
    model = Models(**_gaussian_parameters())
    _, (psi, rho) = model.call_model(
        field_components=1,
        parameters_simulation={
            "gridlength": 4.0, "resol": 4, "lambda_value": 0.0,
            "num_threads": 1, "Boverlap": True, "Plim": 5.6, "rmax": 7.6,
        },
    )
    assert psi.shape == rho.shape == (1, 4, 4, 4)
    assert np.isfinite(psi).all() and np.isfinite(rho).all()


def test_short_gaussian_run_can_be_read_and_visualized(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    output = tmp_path / "gaussian"
    euskera.evolve(
        _gaussian_parameters(),
        evolution_config=euskera.EvolutionConfig(
            gridlength=4.0, resol=4, tmax=0.01
        ),
        output_config=OutputConfig(
            address=str(output),
            save_number=1,
            data_save={
                "grid": True, "save_line_rho": True,
                "save_line_psi": False, "save_line_phi": False,
                "save_energies": False,
            },
        ),
        diagnostics_config={"Numb_Part": False, "Energ": False},
    )

    saved = sorted(output.glob("end_save_line_rho*.npz"))
    assert saved
    with np.load(saved[0]) as archive:
        arrays = [archive[name] for name in archive.files]
    assert arrays and all(np.isfinite(array).all() for array in arrays)
    line = arrays[-1]

    figure, _ = imagshow2D(
        line[None, :], contour=False, xlim=(0, line.size),
        ylim=(0, 1), out=True,
    )
    assert figure.axes
    plt.close(figure)
