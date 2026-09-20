"""Regression checks for explicit active-field boundary corrections."""

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from euskera.backgrounds.fitting import algebSyst, fitting
from euskera.backgrounds.systems import systemMultFreqTot


@pytest.mark.parametrize("active", [None, [True, True]])
def test_zero_rhs_does_not_remove_a_coupled_equation(active):
    matrix = np.array([[1., 1.], [1., 2.]])
    result = algebSyst([matrix, [0., -1.], [0., 0.], [0., 0.]], active=active)
    np.testing.assert_allclose(result, [-1., 1.])
    np.testing.assert_allclose(matrix @ result, [0., 1.])


def test_inactive_field_can_be_between_active_fields(capsys):
    matrix = np.array([[1., 0., 1.], [0., 0., 0.], [1., 0., 2.]])
    result = algebSyst(
        [matrix, [0., 0., -1.], np.zeros(3), np.zeros(3)],
        active=[True, False, True], info=True,
    )
    np.testing.assert_allclose(result, [-1., 0., 1.])
    assert "True" in capsys.readouterr().out


def test_all_inactive_fields_return_zero():
    result = algebSyst(
        [np.zeros((2, 2)), np.zeros(2), np.zeros(2), np.zeros(2)],
        active=[False, False],
    )
    np.testing.assert_array_equal(result, [0., 0.])


@pytest.mark.parametrize("active", [[True], [[True, False]], [1, 0], True])
def test_rejects_invalid_masks(active):
    with pytest.raises(ValueError, match="active"):
        algebSyst([np.eye(2), np.zeros(2), np.zeros(2), np.zeros(2)], active=active)


def test_singular_active_system_is_not_silently_reduced():
    with pytest.raises(np.linalg.LinAlgError):
        algebSyst([np.zeros((2, 2)), np.zeros(2), np.zeros(2), np.zeros(2)])


def test_fitting_preserves_full_background_state_with_inactive_field():
    # Manufacture a reachable boundary from the actual 48-variable ODE, then
    # recover it from perturbed potential seeds. The z field stays zero.
    initial = np.zeros(48)
    initial[[0, 2]] = [0.3, 0.2]
    initial[[6, 8]] = [0.4, 0.6]
    initial[[18, 32, 46]] = 1.0
    limit = [0., 0.4]
    target = solve_ivp(
        systemMultFreqTot, limit, initial, args=([0.],),
        rtol=1e-11, atol=1e-13,
    )
    assert target.success
    boundary = target.y[[0, 2, 4], -1]
    guess = initial.copy()
    guess[[6, 8]] += [0.02, -0.03]
    indck = np.zeros(48, dtype=bool)
    indck[[6, 8, 10]] = True
    indXc = np.zeros(48, dtype=bool)
    indXc[[0, 2, 4]] = True
    inddXc = np.zeros(48, dtype=bool)
    inddXc[[12, 14, 16, 24, 26, 28, 36, 38, 40]] = True

    result = fitting(
        systemMultFreqTot, guess, indck, indXc, boundary, inddXc,
        limit, argf=[0.], active=[True, True, False],
        tol=1e-10, Rtol=1e-11, Atol=1e-13, klim=10,
    )
    assert result.shape == (48,)
    np.testing.assert_array_equal(result[~indck], initial[~indck])
    np.testing.assert_allclose(result[indck], initial[indck], atol=1e-7)
    final = solve_ivp(
        systemMultFreqTot, limit, result, args=([0.],),
        rtol=1e-11, atol=1e-13,
    )
    assert final.success
    np.testing.assert_allclose(final.y[[0, 2, 4], -1], boundary, atol=1e-10)
    np.testing.assert_array_equal(final.y[[4, 5]], 0.)


def test_fitting_validates_mask_even_if_initial_state_already_converged():
    with pytest.raises(ValueError, match="active"):
        fitting(lambda r, y, arg: np.zeros_like(y), [0.], [True], [True],
                [0.], [False], [0., 1.], active=[1])
