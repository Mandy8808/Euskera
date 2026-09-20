"""Numerical regression tests for the standalone spectral package."""

import numpy as np

from euskera.spectral import backgroundOper, cheb, chevQuant
from euskera.spectral.blocks import linBlock


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
