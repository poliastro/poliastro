# coding: utf-8
import pytest

import numpy as np
from numpy.testing import (assert_equal, assert_almost_equal,
                           assert_array_almost_equal)

from astropy import units as u
from astropy import time

from poliastro.bodies import Sun, Earth
from poliastro.twobody import State


def test_state_has_attractor_given_in_constructor():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    ss = State.from_classical(Sun, _d, _, _a, _a, _a, _a)
    assert ss.attractor == Sun


def test_default_time_for_new_state():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    _body = Sun  # Unused body
    expected_epoch = time.Time("J2000", scale='utc')
    ss = State.from_classical(_body, _d, _, _a, _a, _a, _a)
    assert ss.epoch == expected_epoch


def test_state_has_elements_given_in_constructor():
    # Mars data from HORIZONS at J2000
    a = 1.523679 * u.AU
    ecc = 0.093315 * u.one
    p = a * (1 - ecc**2)
    inc = 1.85 * u.deg
    raan = 49.562 * u.deg
    argp = 286.537 * u.deg
    nu = 23.33 * u.deg
    ss = State.from_classical(Sun, p, ecc, inc, raan, argp, nu)
    assert ss.coe() == (p, ecc, inc, raan, argp, nu)


def test_state_has_individual_elements():
    a = 1.523679 * u.AU
    ecc = 0.093315 * u.one
    p = a * (1 - ecc**2)
    inc = 1.85 * u.deg
    raan = 49.562 * u.deg
    argp = 286.537 * u.deg
    nu = 23.33 * u.deg
    ss = State.from_classical(Sun, p, ecc, inc, raan, argp, nu)
    assert ss.a == a
    assert ss.ecc == ecc
    assert ss.inc == inc
    assert ss.raan == raan
    assert ss.argp == argp
    assert ss.nu == nu


def test_state_raises_unitserror_if_elements_units_are_wrong():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    wrong_angle = 1.0 * u.AU
    with pytest.raises(u.UnitsError) as excinfo:
        ss = State.from_classical(Sun, _d, _, _a, _a, _a, wrong_angle)
    assert ("UnitsError: Units must be consistent"
            in excinfo.exconly())


def test_state_has_rv_given_in_constructor():
    r = [1.0, 0.0, 0.0] * u.AU
    v = [0.0, 1.0e-6, 0.0] * u.AU / u.s
    ss = State.from_vectors(Sun, r, v)
    assert ss.rv() == (r, v)


def test_state_raises_unitserror_if_rv_units_are_wrong():
    _d = [1.0, 0.0, 0.0] * u.AU
    wrong_v = [0.0, 1.0e-6, 0.0] * u.AU
    with pytest.raises(u.UnitsError) as excinfo:
        ss = State.from_vectors(Sun, _d, wrong_v)
    assert ("UnitsError: Units must be consistent"
            in excinfo.exconly())


def test_circular_has_proper_semimajor_axis():
    alt = 500 * u.km
    attractor = Earth
    expected_a = Earth.R + alt
    ss = State.circular(attractor, alt)
    assert ss.a == expected_a


def test_geosync_has_proper_period():
    expected_period = 1436  # min
    ss = State.circular(Earth, alt=42164 * u.km - Earth.R)
    assert_almost_equal(ss.period.to(u.min).value, expected_period, decimal=0)


def test_parabolic_has_proper_eccentricity():
    attractor = Earth
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    expected_ecc = 1.0 * u.one
    ss = State.parabolic(attractor, _d, _a, _a, _a, _a)
    assert_almost_equal(ss.ecc, expected_ecc)


def test_parabolic_has_infinite_semimajor_axis():
    attractor = Earth
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    expected_a = np.inf * u.km
    ss = State.parabolic(attractor, _d, _a, _a, _a, _a)
    assert_equal(ss.a, expected_a)


def test_parabolic_has_zero_energy():
    attractor = Earth
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    ss = State.parabolic(attractor, _d, _a, _a, _a, _a)
    assert_almost_equal(ss.energy.value, 0.0)


def test_perigee_and_apogee():
    expected_r_a = 500 * u.km
    expected_r_p = 300 * u.km
    a = (expected_r_a + expected_r_p) / 2
    ecc = expected_r_a / a - 1
    p = a * (1 - ecc**2)
    _a = 1.0 * u.deg  # Unused angle
    ss = State.from_classical(Earth, p, ecc, _a, _a, _a, _a)
    assert_almost_equal(ss.r_a.to(u.km).value,
                        expected_r_a.to(u.km).value)
    assert_almost_equal(ss.r_p.to(u.km).value,
                        expected_r_p.to(u.km).value)


def test_convert_from_rv_to_coe():
    # Data from Vallado, example 2.6
    attractor = Earth
    p = 11067.790 * u.km
    ecc = 0.83285 * u.one
    inc = 87.87 * u.deg
    raan = 227.89 * u.deg
    argp = 53.38 * u.deg
    nu = 92.335 * u.deg
    expected_r = [6525.344, 6861.535, 6449.125]  # km
    expected_v = [4.902276, 5.533124, -1.975709]  # km / s
    r, v = State.from_classical(Earth, p, ecc, inc, raan, argp, nu).rv()
    assert_array_almost_equal(r.value, expected_r, decimal=1)
    assert_array_almost_equal(v.value, expected_v, decimal=5)


def test_convert_from_coe_to_rv():
    # Data from Vallado, example 2.5
    attractor = Earth
    r = [6524.384, 6862.875, 6448.296] * u.km
    v = [4.901327, 5.533756, -1.976341] * u.km / u.s
    ss = State.from_vectors(Earth, r, v)
    p, ecc, inc, raan, argp, nu = ss.coe()
    assert_almost_equal(p.to(u.km).value, 11067.79, decimal=0)
    assert_almost_equal(ecc.value, 0.832853, decimal=4)
    assert_almost_equal(inc.to(u.deg).value, 87.870, decimal=2)
    assert_almost_equal(raan.to(u.deg).value, 227.89, decimal=1)
    assert_almost_equal(argp.to(u.deg).value, 53.38, decimal=2)
    assert_almost_equal(nu.to(u.deg).value, 92.335, decimal=2)


def test_apply_zero_maneuver_returns_equal_state():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    ss = State.from_classical(Sun, _d, _, _a, _a, _a, _a)
    dt = 0 * u.s
    dv = [0, 0, 0] * u.km / u.s
    ss_new = ss.apply_maneuver([(dt, dv)])
    assert_almost_equal(ss_new.r.to(u.km).value,
                        ss.r.to(u.km).value)
    assert_almost_equal(ss_new.v.to(u.km / u.s).value,
                        ss.v.to(u.km / u.s).value)


def test_apply_maneuver_changes_epoch():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    ss = State.from_classical(Sun, _d, _, _a, _a, _a, _a)
    dt = 1 * u.h
    dv = [0, 0, 0] * u.km / u.s
    ss_new = ss.apply_maneuver([(dt, dv)])
    assert ss_new.epoch == ss.epoch + dt


def test_perifocal_points_to_perigee():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    ss = State.from_classical(Sun, _d, _, _a, _a, _a, _a)
    p, _, _ = ss.pqw()
    assert_almost_equal(p, ss.e_vec / ss.ecc)


def test_pqw_for_circular_equatorial_orbit():
    ss = State.circular(Earth, 600 * u.km)
    expected_p = [1, 0, 0] * u.one
    expected_q = [0, 1, 0] * u.one
    expected_w = [0, 0, 1] * u.one
    p, q, w = ss.pqw()
    assert_almost_equal(p, expected_p)
    assert_almost_equal(q, expected_q)
    assert_almost_equal(w, expected_w)
