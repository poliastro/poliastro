from pytest import approx
import pytest

import numpy as np
from numpy.testing import assert_allclose

from astropy import units as u
from astropy import time
from astropy.tests.helper import assert_quantity_allclose

from poliastro.twobody.rv import rv2coe
from poliastro.constants import J2000
from poliastro.bodies import Sun, Earth
from poliastro.twobody import Orbit
from poliastro.twobody.propagation import cowell, kepler, mean_motion
from poliastro.examples import iss

from poliastro.neos import dastcom5

from poliastro.util import norm


@pytest.mark.parametrize('method', [kepler, mean_motion, cowell])
def test_propagation(method):
    # Data from Vallado, example 2.4
    r0 = [1131.340, -2282.343, 6672.423] * u.km
    v0 = [-5.64305, 4.30333, 2.42879] * u.km / u.s
    expected_r = [-4219.7527, 4363.0292, -3958.7666] * u.km
    expected_v = [3.689866, -1.916735, -6.112511] * u.km / u.s

    ss0 = Orbit.from_vectors(Earth, r0, v0)
    tof = 40 * u.min
    ss1 = ss0.propagate(tof, method=method)

    r, v = ss1.rv()

    assert_quantity_allclose(r, expected_r, rtol=1e-5)
    assert_quantity_allclose(v, expected_v, rtol=1e-4)


def test_propagating_to_certain_nu_is_correct():
    # take an elliptic orbit
    a = 1.0 * u.AU
    ecc = 1.0 / 3.0 * u.one
    _a = 0.0 * u.rad
    elliptic = Orbit.from_classical(Sun, a, ecc, _a, _a, _a, _a)
    r_ini, _ = elliptic.rv()

    elliptic_at_perihelion = elliptic.propagate(0.0 * u.rad)
    r_per, _ = elliptic_at_perihelion.rv()

    elliptic_at_aphelion = elliptic.propagate(np.pi * u.rad)
    r_ap, _ = elliptic_at_aphelion.rv()

    assert_quantity_allclose(r_per, r_ini)
    assert_quantity_allclose(norm(r_per), a * (1.0 - ecc))
    assert_quantity_allclose(norm(r_ap), a * (1.0 + ecc))

    # test 10 random true anomaly values
    for _ in range(10):
        nu = np.random.uniform(low=0.0, high=2 * np.pi)
        elliptic = elliptic.propagate(nu * u.rad)
        r, _ = elliptic.rv()
        assert_quantity_allclose(norm(r), a * (1.0 - ecc ** 2) / (1 + ecc * np.cos(nu)))


def test_propagate_accepts_timedelta():
    # Data from Vallado, example 2.4
    r0 = [1131.340, -2282.343, 6672.423] * u.km
    v0 = [-5.64305, 4.30333, 2.42879] * u.km / u.s
    expected_r = [-4219.7527, 4363.0292, -3958.7666] * u.km
    expected_v = [3.689866, -1.916735, -6.112511] * u.km / u.s

    ss0 = Orbit.from_vectors(Earth, r0, v0)
    tof = time.TimeDelta(40 * u.min)
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


def test_propagation_mean_motion_parabolic():
    # example from Howard Curtis (3rd edition), section 3.5, problem 3.15
    p = 2.0 * 6600 * u.km
    _a = 0.0 * u.deg
    orbit = Orbit.parabolic(Earth, p, _a, _a, _a, _a)
    orbit = orbit.propagate(0.8897 / 2.0 * u.h, method=mean_motion)

    _, _, _, _, _, nu0 = rv2coe(Earth.k.to(u.km**3 / u.s**2).value,
                                orbit.r.to(u.km).value,
                                orbit.v.to(u.km / u.s).value)
    assert_quantity_allclose(nu0, np.deg2rad(90.0), rtol=1e-4)

    orbit = Orbit.parabolic(Earth, p, _a, _a, _a, _a)
    orbit = orbit.propagate(36.0 * u.h, method=mean_motion)
    assert_quantity_allclose(norm(orbit.r.to(u.km).value), 304700.0, rtol=1e-4)


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

    r0 = np.array([1131.340, -2282.343, 6672.423])  # km
    v0 = np.array([-5.64305, 4.30333, 2.42879])  # km/s
    tof = 40 * 60.0  # s
    orbit = Orbit.from_vectors(Earth, r0 * u.km, v0 * u.km / u.s)

    results = []

    def cb(t, u_):
        row = [t]
        row.extend(u_)
        results.append(row)

    r, v = cowell(orbit, tof, callback=cb)

    assert len(results) == 17
    assert len(results[0]) == 7
    assert results[-1][0] == tof


def test_cowell_propagation_with_zero_acceleration_equals_kepler():
    # Data from Vallado, example 2.4

    r0 = np.array([1131.340, -2282.343, 6672.423])  # km
    v0 = np.array([-5.64305, 4.30333, 2.42879])  # km/s
    tof = 40 * 60.0  # s

    orbit = Orbit.from_vectors(Earth, r0 * u.km, v0 * u.km / u.s)

    expected_r = np.array([-4219.7527, 4363.0292, -3958.7666])
    expected_v = np.array([3.689866, -1.916735, -6.112511])

    r, v = cowell(orbit, tof, ad=None)

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

    r, v = cowell(ss,
                  tof.to(u.s).value,
                  ad=constant_accel)

    ss_final = Orbit.from_vectors(
        Earth, r * u.km, v * u.km / u.s)

    da_a0 = (ss_final.a - ss.a) / ss.a
    dv_v0 = abs(norm(ss_final.v) - norm(ss.v)) / norm(ss.v)
    assert_quantity_allclose(da_a0, 2 * dv_v0, rtol=1e-2)

    dv = abs(norm(ss_final.v) - norm(ss.v))
    accel_dt = accel * u.km / u.s**2 * tof
    assert_quantity_allclose(dv, accel_dt, rtol=1e-2)


def test_propagate_to_date_has_proper_epoch():
    # Data from Vallado, example 2.4
    r0 = [1131.340, -2282.343, 6672.423] * u.km
    v0 = [-5.64305, 4.30333, 2.42879] * u.km / u.s
    init_epoch = J2000
    final_epoch = time.Time("2000-01-01 12:40:00", scale="tdb")

    expected_r = [-4219.7527, 4363.0292, -3958.7666] * u.km
    expected_v = [3.689866, -1.916735, -6.112511] * u.km / u.s

    ss0 = Orbit.from_vectors(Earth, r0, v0, epoch=init_epoch)
    ss1 = ss0.propagate(final_epoch)

    r, v = ss1.rv()

    assert_quantity_allclose(r, expected_r, rtol=1e-5)
    assert_quantity_allclose(v, expected_v, rtol=1e-4)

    # Tolerance should be higher, see https://github.com/astropy/astropy/issues/6638
    assert (ss1.epoch - final_epoch).sec == approx(0.0, abs=1e-6)


@pytest.mark.filterwarnings('ignore:ERFA')
@pytest.mark.filterwarnings('ignore::UserWarning')
@pytest.mark.parametrize('method', [
    mean_motion,
    pytest.param(cowell, marks=pytest.mark.xfail),
    pytest.param(kepler, marks=pytest.mark.xfail),
])
def test_propagate_long_times_keeps_geometry(method):
    # See https://github.com/poliastro/poliastro/issues/265
    time_of_flight = 100 * u.year

    res = iss.propagate(time_of_flight, method=method)

    assert_quantity_allclose(iss.a, res.a)
    assert_quantity_allclose(iss.ecc, res.ecc)
    assert_quantity_allclose(iss.inc, res.inc)
    assert_quantity_allclose(iss.raan, res.raan)
    assert_quantity_allclose(iss.argp, res.argp)

    assert_quantity_allclose((res.epoch - iss.epoch).to(time_of_flight.unit), time_of_flight)


def test_long_propagations_kepler_agrees_mean_motion():
    tof = 100 * u.year
    r_mm, v_mm = iss.propagate(tof, method=mean_motion).rv()
    r_k, v_k = iss.propagate(tof, method=kepler).rv()
    assert_quantity_allclose(r_mm, r_k, rtol=1e-6)
    assert_quantity_allclose(v_mm, v_k, rtol=1e-6)

    halleys = dastcom5.orbit_from_name('1P')[0]
    r_mm, v_mm = halleys.propagate(tof, method=mean_motion).rv()
    r_k, v_k = halleys.propagate(tof, method=kepler).rv()
    assert_quantity_allclose(r_mm, r_k, rtol=1e-6)
    assert_quantity_allclose(v_mm, v_k, rtol=1e-6)
