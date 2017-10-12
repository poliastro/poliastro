import pytest

import numpy as np
from numpy.testing import assert_allclose

from astropy import units as u
from astropy.time import Time
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
    assert_allclose(result, expected_vector, atol=1e-16)


def test_rotate_vector_with_units():
    vector = [1, 0, 0] * u.m
    angle = 90 * u.deg
    axis = 'y'
    expected_vector = [0, 0, -1] * u.m
    result = util.rotate(vector, angle, axis)
    assert_quantity_allclose(result, expected_vector, atol=1e-16 * u.m)


def test_transform_unitless_vector():
    vector = [0, 1, 0]
    angle = 45 * u.deg
    axis = 'z'
    expected_vector = [np.sqrt(2) / 2, np.sqrt(2) / 2, 0]
    result = util.transform(vector, angle, axis)
    assert_allclose(result, expected_vector)


def test_time_range_spacing_periods():
    start_time = "2017-10-12 00:00:00"
    end_time = "2017-10-12 00:04:00"
    spacing = 1 * u.minute
    periods = 5

    expected_scale = "utc"
    expected_duration = 4 * u.min

    result_1 = util.time_range(start_time, spacing=spacing, periods=periods)
    result_2 = util.time_range(start_time, end=end_time, periods=periods)
    result_3 = util.time_range(Time(start_time), end=Time(end_time), periods=periods)

    assert len(result_1) == len(result_2) == len(result_3) == periods
    assert result_1.scale == result_2.scale == result_3.scale == expected_scale

    assert_quantity_allclose((result_1[-1] - result_1[0]).to(u.s), expected_duration)
    assert_quantity_allclose((result_2[-1] - result_2[0]).to(u.s), expected_duration)
    assert_quantity_allclose((result_3[-1] - result_3[0]).to(u.s), expected_duration)


def test_time_range_requires_keyword_arguments():
    with pytest.raises(TypeError) as excinfo:
        util.time_range(0, 0)
    assert "TypeError: time_range() takes 1 positional argument but" in excinfo.exconly()


def test_time_range_raises_error_wrong_arguments():
    exception_message = "ValueError: Either 'end' or 'spacing' must be specified"

    with pytest.raises(ValueError) as excinfo_1:
        util.time_range("2017-10-12 00:00")

    with pytest.raises(ValueError) as excinfo_2:
        util.time_range("2017-10-12 00:00", spacing=0, end=0, periods=0)

    assert exception_message in excinfo_1.exconly()
    assert exception_message in excinfo_2.exconly()
