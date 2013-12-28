import numpy as np
from numpy.testing import assert_array_almost_equal

from astropy import units as u

from poliastro import util


def test_rotate_unitless_vector():
    vector = [1, 0, 0]
    angle = 90 * u.deg
    axis = 'z'
    expected_vector = [0, 1, 0]
    result = util.rotate(vector, angle, axis)
    assert_array_almost_equal(result, expected_vector)


def test_rotate_vector_with_units():
    vector = [1, 0, 0] * u.m
    angle = 90 * u.deg
    axis = 'y'
    expected_vector = [0, 0, -1] * u.m
    result = util.rotate(vector, angle, axis)
    assert_array_almost_equal(result, expected_vector)


def test_transform_unitless_vector():
    vector = [0, 1, 0]
    angle = 45 * u.deg
    axis = 'z'
    expected_vector = [np.sqrt(2) / 2, np.sqrt(2) / 2, 0]
    result = util.transform(vector, angle, axis)
    assert_array_almost_equal(result, expected_vector)
