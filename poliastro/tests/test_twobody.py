# coding: utf-8
import pytest

from numpy.testing import assert_almost_equal, assert_array_almost_equal

from astropy import units as u
from astropy import time

from poliastro.bodies import Sun, Earth
from poliastro import twobody


def test_state_raises_valueerror_if_meaningless_state():
    _ = 1.0
    with pytest.raises(ValueError) as excinfo:
        ss = twobody.State(Sun, (_,))
    assert ("ValueError: Incorrect number of parameters"
            in excinfo.exconly())


def test_state_has_attractor_given_in_constructor():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    ss = twobody.State(Sun, (_d, _, _a, _a, _a, _a))
    assert ss.attractor == Sun


def test_default_time_for_new_state():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    _body = Sun  # Unused body
    expected_epoch = time.Time("J2000", scale='utc')
    ss = twobody.State(_body, (_d, _, _a, _a, _a, _a))
    assert ss.epoch == expected_epoch


def test_state_has_elements_given_in_constructor():
    # Mars data from HORIZONS at J2000
    a = 1.523679 * u.AU
    ecc = 0.093315
    inc = 1.85 * u.deg
    raan = 49.562 * u.deg
    argp = 286.537 * u.deg
    nu = 23.33 * u.deg
    ss = twobody.State(Sun, (a, ecc, inc, raan, argp, nu))
    assert ss.elements == (a, ecc, inc, raan, argp, nu)


def test_state_raises_valueerror_if_elements_units_are_missing():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    wrong_angle = 1.0
    with pytest.raises(ValueError) as excinfo:
        ss = twobody.State(Sun, (_d, _, _a, _a, _a, wrong_angle))
    assert ("Elements must have units (use astropy.units)"
            in excinfo.exconly())


def test_state_raises_unitserror_if_elements_units_are_wrong():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    wrong_angle = 1.0 * u.AU
    with pytest.raises(u.UnitsError) as excinfo:
        ss = twobody.State(Sun, (_d, _, _a, _a, _a, wrong_angle))
    assert ("UnitsError: Units must be consistent"
            in excinfo.exconly())


def test_state_has_rv_given_in_constructor():
    r = [1.0, 0.0, 0.0] * u.AU
    v = [0.0, 1.0e-6, 0.0] * u.AU / u.s
    ss = twobody.State(Sun, (r, v))
    assert ss.rv == (r, v)


def test_state_raises_valueerror_if_rv_units_are_missing():
    _d = [1.0, 0.0, 0.0] * u.AU
    wrong_v = [0.0, 1.0e-6, 0.0]
    with pytest.raises(ValueError) as excinfo:
        ss = twobody.State(Sun, (_d, wrong_v))
    assert ("ValueError: r and v vectors must have units (use astropy.units)"
            in excinfo.exconly())


def test_state_raises_unitserror_if_rv_units_are_wrong():
    _d = [1.0, 0.0, 0.0] * u.AU
    wrong_v = [0.0, 1.0e-6, 0.0] * u.AU
    with pytest.raises(u.UnitsError) as excinfo:
        ss = twobody.State(Sun, (_d, wrong_v))
    assert ("UnitsError: Units must be consistent"
            in excinfo.exconly())


def test_convert_from_rv_to_coe():
    # Data from Vallado, example 2.6
    attractor = Earth
    p = 11067.790 * u.km
    ecc = 0.83285
    a = p / (1 - ecc ** 2)
    inc = 87.87 * u.deg
    raan = 227.89 * u.deg
    argp = 53.38 * u.deg
    nu = 92.335 * u.deg
    expected_r = [6525.344, 6861.535, 6449.125] * u.km
    expected_v = [4.902276, 5.533124, -1.975709] * u.km / u.s
    r, v = twobody.State(Earth, (a, ecc, inc, raan, argp, nu)).rv
    assert_array_almost_equal(r, expected_r, decimal=1)
    assert_array_almost_equal(v, expected_v, decimal=3)


def test_convert_from_coe_to_rv():
    # Data from Vallado, example 2.5
    attractor = Earth
    r = [6524.384, 6862.875, 6448.296] * u.km
    v = [4.901327, 5.533756, -1.976341] * u.km / u.s
    a, ecc, inc, omega, argp, nu = twobody.State(Earth, (r, v)).elements
    p = a * (1 - ecc ** 2)
    assert_almost_equal(p, 11067.79 * u.km)
    assert_almost_equal(ecc, 0.832853, decimal=3)
    assert_almost_equal(inc, 87.870 * u.deg, decimal=2)
    assert_almost_equal(omega, 227.89 * u.deg, decimal=2)
    assert_almost_equal(argp, 53.38 * u.deg, decimal=2)
    assert_almost_equal(nu, 92.335 * u.deg, decimal=1)
