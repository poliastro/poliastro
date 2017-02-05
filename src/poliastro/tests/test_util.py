# coding: utf-8
import numpy as np
from numpy.testing import assert_array_almost_equal

from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from poliastro import util


def test_simple_circular_velocity():
    k = 398600 * u.km ** 3 / u.s ** 2
    a = 7000 * u.km

    expected_V = 7.5460491 * u.km / u.s

    V = util.circular_velocity(k, a)

    assert_quantity_allclose(V, expected_V)


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
    assert_array_almost_equal(result.value, expected_vector.value)


def test_transform_unitless_vector():
    vector = [0, 1, 0]
    angle = 45 * u.deg
    axis = 'z'
    expected_vector = [np.sqrt(2) / 2, np.sqrt(2) / 2, 0]
    result = util.transform(vector, angle, axis)
    assert_array_almost_equal(result, expected_vector)
