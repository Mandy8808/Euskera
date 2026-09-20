"""Tests for shared simulation utilities."""

import pytest

from euskera.observables.config import DiagnosticsConfig
from euskera.observables.simulation_conserv_quant import Conserv
from euskera.tools.tools import dtime, update_simulation_parameters


def test_ji_diagnostics_are_rejected_by_typed_configuration():
    with pytest.raises(ValueError, match="Ji diagnostics are not implemented"):
        DiagnosticsConfig(Ji=True)


def test_ji_diagnostics_are_rejected_by_legacy_configuration():
    with pytest.raises(ValueError, match="Ji diagnostics are not implemented"):
        update_simulation_parameters({"Ji": True}, {"Ji": False})


def test_ji_diagnostics_are_rejected_by_conserv():
    with pytest.raises(ValueError, match="Ji diagnostics are not implemented"):
        Conserv(None, {"Ji": True}, {})


def test_dtime_does_not_add_a_save_interval_when_already_divisible():
    _, iterations_per_save, actual_steps = dtime(
        tmax=3.0,
        gridlength=1.0,
        resol=1,
        step_factor=1.0,
        save_number=5,
    )

    assert actual_steps == 10
    assert iterations_per_save == 2


def test_dtime_rounds_up_to_the_next_save_interval_when_needed():
    _, iterations_per_save, actual_steps = dtime(
        tmax=2.0,
        gridlength=1.0,
        resol=1,
        step_factor=1.0,
        save_number=5,
    )

    assert actual_steps == 10
    assert iterations_per_save == 2
