"""Numerical regression tests for the standalone spectral package."""

import numpy as np

from euskera.spectral import backgroundOper, cheb, chevQuant
from euskera.spectral.blocks import linBlock
from euskera.spectral.eigensolver import _eigenvalue_residuals, spectrum, LamJval


def test_cheb_differentiates_a_quadratic():
    derivative, points = cheb(8)

    np.testing.assert_allclose(derivative @ (points**2), 2 * points, atol=1e-12)


def test_chev_quant_maps_and_scales_the_radial_grid():
    radii, derivative, points = chevQuant([0, 6, 3.0])

    np.testing.assert_allclose(radii[[0, -1]], [0.0, 3.0])
    np.testing.assert_allclose(radii, (1.0 - points) * 1.5)
    np.testing.assert_allclose(
        derivative @ radii, -1.5 * np.ones_like(radii), atol=1e-10
    )


def test_background_operators_have_consistent_dimensions():
    sigma = [lambda radius: np.ones_like(radius)]
    potential = [lambda radius: np.zeros_like(radius)]

    sigma_matrices, potential_matrices, inverse, radial, second, radii, lam_star = (
        backgroundOper((sigma, potential), [1, 6, 4.0], (2.0, -1.0), 0.5)
    )

    assert sigma_matrices[0].shape == (5, 5)
    assert potential_matrices[0].shape == (5, 5)
    assert inverse.shape == (5, 5)
    assert radial.shape == (5, 5)
    assert second.shape == (5, 5)
    assert radii.shape == (5,)
    assert lam_star == 1.75

    np.testing.assert_allclose(
        (second - 2 * radial) @ inverse, np.eye(5), atol=1e-10
    )


def test_linear_block_places_the_expected_operator_terms():
    size = 3
    sigma = np.diag([1.0, 2.0, 3.0])
    potential = np.diag([4.0, 5.0, 6.0])
    inverse = np.eye(size)
    second = 2.0 * np.eye(size)

    blocks = linBlock(
        size + 1,
        [sigma],
        [potential],
        inverse,
        second,
        (2.0, 3.0),
        2.0,
    )
    m11, m12, m21, m22, m33 = blocks

    np.testing.assert_allclose(m12, 0.0)
    np.testing.assert_allclose(m21, 0.0)
    np.testing.assert_allclose(m33, m22)
    np.testing.assert_allclose(
        m11[:size, size:],
        second + potential - sigma @ sigma,
    )
    np.testing.assert_allclose(
        m11[size:, :size],
        second + potential - 3 * sigma @ sigma - 2 * sigma @ inverse @ sigma,
    )


def test_eigenvalue_residuals_are_near_zero_for_exact_eigenpairs():
    operator = np.diag([1.0, 2.0, 3.0]).astype(complex)
    eigenvalues = np.array([1.0, 2.0, 3.0], dtype=complex)
    eigenvectors = np.eye(3, dtype=complex)

    diagnostics = _eigenvalue_residuals(operator, eigenvalues, eigenvectors)

    assert np.all(diagnostics["allclose"])
    np.testing.assert_allclose(diagnostics["absolute_residuals"], 0.0, atol=1e-12)
    np.testing.assert_allclose(diagnostics["relative_residuals"], 0.0, atol=1e-12)
    np.testing.assert_allclose(diagnostics["scaled_residuals"], 0.0, atol=1e-12)


def test_eigenvalue_residuals_flag_an_incorrect_eigenpair():
    operator = np.diag([1.0, 2.0, 3.0]).astype(complex)
    eigenvalues = np.array([1.0, 2.0, 3.0], dtype=complex)
    eigenvectors = np.eye(3, dtype=complex)
    eigenvectors[:, 1] = [0.0, 1.0, 1.0]  # not an eigenvector of the operator

    diagnostics = _eigenvalue_residuals(operator, eigenvalues, eigenvectors)

    assert diagnostics["allclose"][0]
    assert not diagnostics["allclose"][1]
    assert diagnostics["allclose"][2]
    assert diagnostics["absolute_residuals"][1] > 1e-3


def test_eigenvalue_residuals_track_the_size_of_the_perturbation():
    operator = np.diag([1.0, 2.0, 3.0]).astype(complex)
    eigenvalues = np.array([1.0, 2.0, 3.0], dtype=complex)

    tiny_perturbation = np.eye(3, dtype=complex)
    tiny_perturbation[0, 1] = 1e-10

    large_perturbation = np.eye(3, dtype=complex)
    large_perturbation[0, 1] = 1e-2

    tiny_diagnostics = _eigenvalue_residuals(operator, eigenvalues, tiny_perturbation)
    large_diagnostics = _eigenvalue_residuals(operator, eigenvalues, large_perturbation)

    # A perturbation below the default tolerance still passes the boolean check...
    assert tiny_diagnostics["allclose"][1]
    # ...while a much larger one is correctly flagged as inconsistent.
    assert not large_diagnostics["allclose"][1]
    # The quantitative residual grows with the size of the perturbation.
    assert (
        tiny_diagnostics["absolute_residuals"][1]
        < large_diagnostics["absolute_residuals"][1]
    )


def test_lamjval_diagnostics_match_the_reported_spectrum():
    functions = ([lambda r: np.exp(-r * r)], [lambda r: -0.3 * np.ones_like(r)])

    for real_only in (False, True):
        spectrum, extra, diagnostics = LamJval(
            functions, 4.0, (1.0, 0.2), 0, "linear",
            Nptos=6, Jval=[0, 1], real=real_only, diagnostics=True,
        )

        assert len(diagnostics) == len(spectrum) == len(extra)
        for (j_spectrum, values), (j_diag, residuals) in zip(spectrum, diagnostics):
            assert j_spectrum == j_diag
            assert residuals["allclose"].shape == values.shape
            assert residuals["absolute_residuals"].shape == values.shape
            # The reported eigenpairs should satisfy the eigenproblem well.
            assert np.all(residuals["absolute_residuals"] < 1e-8)


def test_spectrum_without_diagnostics_keeps_the_four_element_return():
    functions = ([lambda r: np.exp(-r * r)], [lambda r: -0.3 * np.ones_like(r)])

    result = spectrum(functions, [0, 6, 4.0], (1.0, 0.2), 0, "linear")

    assert len(result) == 4
