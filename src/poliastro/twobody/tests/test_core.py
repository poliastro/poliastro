# coding: utf-8
import pytest

import numpy as np

from numpy.testing import (assert_almost_equal,
                           assert_array_almost_equal)

from astropy import units as u

from poliastro.bodies import Sun, Earth
from poliastro.twobody import State

from poliastro.twobody.propagation import cowell

from poliastro.util import norm


def test_state_has_attractor_given_in_constructor():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    ss = State.from_classical(Sun, _d, _, _a, _a, _a, _a)
    assert ss.attractor == Sun


def test_state_has_elements_given_in_constructor():
    # Mars data from HORIZONS at J2000
    a = 1.523679 * u.AU
    ecc = 0.093315 * u.one
    inc = 1.85 * u.deg
    raan = 49.562 * u.deg
    argp = 286.537 * u.deg
    nu = 23.33 * u.deg
    ss = State.from_classical(Sun, a, ecc, inc, raan, argp, nu)
    assert ss.coe() == (a, ecc, inc, raan, argp, nu)


def test_state_has_individual_elements():
    a = 1.523679 * u.AU
    ecc = 0.093315 * u.one
    inc = 1.85 * u.deg
    raan = 49.562 * u.deg
    argp = 286.537 * u.deg
    nu = 23.33 * u.deg
    ss = State.from_classical(Sun, a, ecc, inc, raan, argp, nu)
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
    assert ("UnitsError: Argument 'nu' to function 'from_classical' must be in units convertible to 'rad'."
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
    assert ("UnitsError: Argument 'v' to function 'from_vectors' must be in units convertible to 'm / s'."
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


def test_parabolic_elements_fail_early():
    attractor = Earth
    ecc = 1.0 * u.one
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    with pytest.raises(ValueError) as excinfo:
        ss = State.from_classical(attractor, _d, ecc, _a, _a, _a, _a)
    assert ("ValueError: For parabolic orbits use State.parabolic instead"
            in excinfo.exconly())


def test_parabolic_has_proper_eccentricity():
    attractor = Earth
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    expected_ecc = 1.0 * u.one
    ss = State.parabolic(attractor, _d, _a, _a, _a, _a)
    assert_almost_equal(ss.ecc, expected_ecc)


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
    _a = 1.0 * u.deg  # Unused angle
    ss = State.from_classical(Earth, a, ecc, _a, _a, _a, _a)
    assert_almost_equal(ss.r_a.to(u.km).value,
                        expected_r_a.to(u.km).value)
    assert_almost_equal(ss.r_p.to(u.km).value,
                        expected_r_p.to(u.km).value)


def test_convert_from_rv_to_coe():
    # Data from Vallado, example 2.6
    attractor = Earth
    p = 11067.790 * u.km
    ecc = 0.83285 * u.one
    a = p / (1 - ecc ** 2)
    inc = 87.87 * u.deg
    raan = 227.89 * u.deg
    argp = 53.38 * u.deg
    nu = 92.335 * u.deg
    expected_r = [6525.344, 6861.535, 6449.125]  # km
    expected_v = [4.902276, 5.533124, -1.975709]  # km / s
    r, v = State.from_classical(Earth, a, ecc, inc, raan, argp, nu).rv()
    assert_array_almost_equal(r.value, expected_r, decimal=1)
    assert_array_almost_equal(v.value, expected_v, decimal=5)


def test_convert_from_coe_to_rv():
    # Data from Vallado, example 2.5
    attractor = Earth
    r = [6524.384, 6862.875, 6448.296] * u.km
    v = [4.901327, 5.533756, -1.976341] * u.km / u.s
    ss = State.from_vectors(Earth, r, v)
    _, ecc, inc, raan, argp, nu = ss.coe()
    p = ss.p
    assert_almost_equal(p.to(u.km).value, 11067.79, decimal=0)
    assert_almost_equal(ecc.value, 0.832853, decimal=4)
    assert_almost_equal(inc.to(u.deg).value, 87.870, decimal=2)
    assert_almost_equal(raan.to(u.deg).value, 227.89, decimal=1)
    assert_almost_equal(argp.to(u.deg).value, 53.38, decimal=2)
    assert_almost_equal(nu.to(u.deg).value, 92.335, decimal=2)


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


def test_pqw_returns_dimensionless():
    r_0 = ([1, 0, 0] * u.au).to(u.km)
    v_0 = ([0, 6, 0] * u.au / u.year).to(u.km / u.day)
    ss = State.from_vectors(Sun, r_0, v_0)

    p, q, w = ss.pqw()

    assert p.unit == u.one
    assert q.unit == u.one
    assert w.unit == u.one


def test_cowell_propagation_with_zero_acceleration_equals_kepler():
    # Data from Vallado, example 2.4
    k = Earth.k.to(u.km**3 / u.s**2).value

    r0 = np.array([1131.340, -2282.343, 6672.423])  # km
    v0 = np.array([-5.64305, 4.30333, 2.42879])  # km/s
    tof = 40 * 60.0  # s

    expected_r = np.array([-4219.7527, 4363.0292, -3958.7666])
    expected_v = np.array([3.689866, -1.916735, -6.112511])

    r, v = cowell(k, r0, v0, tof, None)

    assert_array_almost_equal(r, expected_r, decimal=1)
    assert_array_almost_equal(v, expected_v, decimal=4)


def test_cowell_propagation_circle_to_circle():
    # From [Edelbaum, 1961]
    accel = 1e-7

    def constant_accel(t0, u, k):
        v = u[3:]
        norm_v = (v[0]**2 + v[1]**2 + v[2]**2)**.5
        return accel * v / norm_v

    ss = State.circular(Earth, 500 * u.km)
    tof = 20 * ss.period

    r0, v0 = ss.rv()
    k = ss.attractor.k

    r, v = cowell(k.to(u.km**3 / u.s**2).value,
                  r0.to(u.km).value,
                  v0.to(u.km / u.s).value,
                  tof.to(u.s).value,
                  constant_accel)

    ss_final = State.from_vectors(Earth,
                                  r * u.km,
                                  v * u.km / u.s)

    da_a0 = (ss_final.a - ss.a) / ss.a
    dv_v0 = abs(norm(ss_final.v) - norm(ss.v)) / norm(ss.v)
    assert_almost_equal(da_a0.value, 2 * dv_v0.value, decimal=4)

    dv = abs(norm(ss_final.v) - norm(ss.v))
    accel_dt = accel * u.km / u.s**2 * tof
    assert_almost_equal(dv.value, accel_dt.value, decimal=4)
