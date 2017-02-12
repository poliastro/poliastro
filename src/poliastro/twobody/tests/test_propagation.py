# coding: utf-8
from numpy.testing import assert_allclose
import numpy as np

from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from poliastro.bodies import Sun, Earth
from poliastro.twobody import Orbit
from poliastro.twobody.propagation import cowell

from poliastro.util import norm


def test_propagation():
    # Data from Vallado, example 2.4
    r0 = [1131.340, -2282.343, 6672.423] * u.km
    v0 = [-5.64305, 4.30333, 2.42879] * u.km / u.s
    expected_r = [-4219.7527, 4363.0292, -3958.7666] * u.km
    expected_v = [3.689866, -1.916735, -6.112511] * u.km / u.s

    ss0 = Orbit.from_vectors(Earth, r0, v0)
    tof = 40 * u.min
    ss1 = ss0.propagate(tof)

    r, v = ss1.rv()

    assert_quantity_allclose(r, expected_r, rtol=1e-5)
    assert_quantity_allclose(v, expected_v, rtol=1e-4)


def test_propagation_hyperbolic():
    # Data from Curtis, example 3.5
    r0 = [Earth.R.to(u.km).value + 300, 0, 0] * u.km
    v0 = [0, 15, 0] * u.km / u.s
    expected_r_norm = 163180 * u.km
    expected_v_norm = 10.51 * u.km / u.s

    ss0 = Orbit.from_vectors(Earth, r0, v0)
    tof = 14941 * u.s
    ss1 = ss0.propagate(tof)
    r, v = ss1.rv()

    assert_quantity_allclose(norm(r), expected_r_norm, rtol=1e-4)
    assert_quantity_allclose(norm(v), expected_v_norm, rtol=1e-3)


def test_propagation_zero_time_returns_same_state():
    # Bug #50
    r0 = [1131.340, -2282.343, 6672.423] * u.km
    v0 = [-5.64305, 4.30333, 2.42879] * u.km / u.s
    ss0 = Orbit.from_vectors(Earth, r0, v0)
    tof = 0 * u.s

    ss1 = ss0.propagate(tof)

    r, v = ss1.rv()

    assert_allclose(r.value, r0.value)
    assert_allclose(v.value, v0.value)


def test_apply_zero_maneuver_returns_equal_state():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    ss = Orbit.from_classical(Sun, _d, _, _a, _a, _a, _a)
    dt = 0 * u.s
    dv = [0, 0, 0] * u.km / u.s
    orbit_new = ss.apply_maneuver([(dt, dv)])
    assert_allclose(orbit_new.r.to(u.km).value,
                        ss.r.to(u.km).value)
    assert_allclose(orbit_new.v.to(u.km / u.s).value,
                        ss.v.to(u.km / u.s).value)


def test_cowell_propagation_callback():
    # Data from Vallado, example 2.4
    k = Earth.k.to(u.km**3 / u.s**2).value

    r0 = np.array([1131.340, -2282.343, 6672.423])  # km
    v0 = np.array([-5.64305, 4.30333, 2.42879])  # km/s
    tof = 40 * 60.0  # s

    results = []

    def cb(t, u_):
        row = [t]
        row.extend(u_)
        results.append(row)

    r, v = cowell(k, r0, v0, tof, callback=cb)

    assert len(results) == 17
    assert len(results[0]) == 7
    assert results[-1][0] == tof


def test_cowell_propagation_with_zero_acceleration_equals_kepler():
    # Data from Vallado, example 2.4
    k = Earth.k.to(u.km**3 / u.s**2).value

    r0 = np.array([1131.340, -2282.343, 6672.423])  # km
    v0 = np.array([-5.64305, 4.30333, 2.42879])  # km/s
    tof = 40 * 60.0  # s

    expected_r = np.array([-4219.7527, 4363.0292, -3958.7666])
    expected_v = np.array([3.689866, -1.916735, -6.112511])

    r, v = cowell(k, r0, v0, tof, ad=None)

    assert_allclose(r, expected_r, rtol=1e-5)
    assert_allclose(v, expected_v, rtol=1e-4)


def test_cowell_propagation_circle_to_circle():
    # From [Edelbaum, 1961]
    accel = 1e-7

    def constant_accel(t0, u, k):
        v = u[3:]
        norm_v = (v[0]**2 + v[1]**2 + v[2]**2)**.5
        return accel * v / norm_v

    ss = Orbit.circular(Earth, 500 * u.km)
    tof = 20 * ss.period

    r0, v0 = ss.rv()
    k = ss.attractor.k

    r, v = cowell(k.to(u.km**3 / u.s**2).value,
                  r0.to(u.km).value,
                  v0.to(u.km / u.s).value,
                  tof.to(u.s).value,
                  ad=constant_accel)

    ss_final = Orbit.from_vectors(Earth,
                       r * u.km,
                       v * u.km / u.s)

    da_a0 = (ss_final.a - ss.a) / ss.a
    dv_v0 = abs(norm(ss_final.v) - norm(ss.v)) / norm(ss.v)
    assert_quantity_allclose(da_a0, 2 * dv_v0, rtol=1e-2)

    dv = abs(norm(ss_final.v) - norm(ss.v))
    accel_dt = accel * u.km / u.s**2 * tof
    assert_quantity_allclose(dv, accel_dt, rtol=1e-2)
