from functools import partial

import numpy as np
import pytest
from astropy.coordinates.matrix_utilities import (
    rotation_matrix as rotation_matrix_astropy,
)
from hypothesis import given, settings, strategies as st
from numpy.testing import assert_allclose, assert_array_equal

from poliastro.core.util import (
    alinspace,
    rotation_matrix as rotation_matrix_poliastro,
    spherical_to_cartesian,
)


def _test_rotation_matrix_with_v(v, angle, axis):
    exp = rotation_matrix_astropy(np.degrees(-angle), "xyz"[axis]) @ v
    res = rotation_matrix_poliastro(angle, axis) @ v
    assert_allclose(exp, res)


def _test_rotation_matrix(angle, axis):
    expected = rotation_matrix_astropy(-np.rad2deg(angle), "xyz"[axis])
    result = rotation_matrix_poliastro(angle, axis)
    assert_allclose(expected, result)


def test_rotation_matrix():
    v = np.array([-0.30387748, -1.4202498, 0.24305655])
    for angle in (0.5, np.array([-np.pi, np.pi])):
        for axis in (0, 1, 2):
            _test_rotation_matrix_with_v(v, angle, axis)


# These tests are adapted from astropy:
# https://github.com/astropy/astropy/blob/main/astropy/coordinates/tests/test_matrix_utilities.py
def test_rotation_matrix_astropy():
    assert_array_equal(rotation_matrix_poliastro(0, 0), np.eye(3))
    assert_allclose(
        rotation_matrix_poliastro(np.deg2rad(-90), 1),
        [[0, 0, -1], [0, 1, 0], [1, 0, 0]],
        atol=1e-12,
    )

    assert_allclose(
        rotation_matrix_poliastro(np.deg2rad(90), 2),
        [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
        atol=1e-12,
    )

    # make sure it also works for very small angles
    assert_allclose(
        rotation_matrix_astropy(-0.000001, "x"),
        rotation_matrix_poliastro(np.deg2rad(0.000001), 0),
    )


def test_rotation_matrix_x():
    _test_rotation_matrix(0.218, 0)


def test_rotation_matrix_y():
    _test_rotation_matrix(0.218, 1)


def test_rotation_matrix_z():
    _test_rotation_matrix(0.218, 2)


def test_spherical_to_cartesian():
    result = spherical_to_cartesian(np.array([0.5, np.pi / 4, -np.pi / 4]))
    expected = np.array([0.25, -0.25, 0.35355339])
    assert np.allclose(expected, result)

    result = spherical_to_cartesian(np.array([0.5, -np.pi / 4, np.pi / 4]))
    expected = np.array([-0.25, -0.25, 0.35355339])
    assert np.allclose(expected, result)

    result = spherical_to_cartesian(
        np.array([[0.5, np.pi / 4, -np.pi / 4], [0.5, -np.pi / 4, np.pi / 4]])
    )
    expected = np.array([[0.25, -0.25, 0.35355339], [-0.25, -0.25, 0.35355339]])
    assert np.allclose(expected, result)


angles = partial(st.floats, min_value=-2 * np.pi, max_value=2 * np.pi)


@settings(deadline=None)
@given(
    x=angles(),
    y=st.one_of(angles(), st.none()),
)
def test_alinspace_is_always_increasing(x, y):
    result = alinspace(x, y)

    assert (np.diff(result) >= 0).all()


@settings(deadline=None)
@given(
    x=st.floats(min_value=-np.pi, max_value=np.pi),
    y=st.one_of(st.floats(min_value=-np.pi, max_value=np.pi), st.none()),
)
def test_alinspace_is_always_increasing_with_angles_inside_anomaly_range(x, y):
    result = alinspace(x, y)

    assert (np.diff(result) >= 0).all()


@settings(deadline=None)
@given(x=angles())
def test_alinspace_no_max_value_uses_full_circle(x):
    result = alinspace(x)

    assert result.max() - result.min() == pytest.approx(2 * np.pi)
