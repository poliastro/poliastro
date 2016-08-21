# coding: utf-8
from numpy.testing import assert_almost_equal, assert_array_almost_equal

from astropy import units as u
from astropy import time

from poliastro.bodies import Sun, Earth
from poliastro.twobody import State
from poliastro.twobody.orbit import Orbit

from poliastro.util import norm


def test_default_time_for_new_state():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    _body = Sun  # Unused body
    expected_epoch = time.Time("J2000", scale='utc')
    ss = State.from_classical(_body, _d, _, _a, _a, _a, _a)
    orbit = Orbit(ss)
    assert orbit.epoch == expected_epoch


def test_apply_maneuver_changes_epoch():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    ss = State.from_classical(Sun, _d, _, _a, _a, _a, _a)
    dt = 1 * u.h
    dv = [0, 0, 0] * u.km / u.s
    orbit = Orbit(ss)
    orbit_new = orbit.apply_maneuver([(dt, dv)])
    assert orbit_new.epoch == orbit.epoch + dt


def test_propagation():
    # Data from Vallado, example 2.4
    r0 = [1131.340, -2282.343, 6672.423] * u.km
    v0 = [-5.64305, 4.30333, 2.42879] * u.km / u.s
    ss0 = State.from_vectors(Earth, r0, v0)
    tof = 40 * u.min
    or0 = Orbit(ss0)
    or1 = or0.propagate(tof)
    r, v = or1.state.rv()
    assert_array_almost_equal(r.value, [-4219.7527, 4363.0292, -3958.7666],
                              decimal=1)
    assert_array_almost_equal(v.value, [3.689866, -1.916735, -6.112511],
                              decimal=4)


def test_propagation_hyperbolic():
    # Data from Curtis, example 3.5
    r0 = [Earth.R.to(u.km).value + 300, 0, 0] * u.km
    v0 = [0, 15, 0] * u.km / u.s
    ss0 = State.from_vectors(Earth, r0, v0)
    tof = 14941 * u.s
    or0 = Orbit(ss0)
    or1 = or0.propagate(tof)
    r, v = or1.state.rv()
    assert_almost_equal(norm(r).to(u.km).value, 163180, decimal=-1)
    assert_almost_equal(norm(v).to(u.km/u.s).value, 10.51, decimal=2)


def test_propagation_zero_time_returns_same_state():
    # Bug #50
    r0 = [1131.340, -2282.343, 6672.423] * u.km
    v0 = [-5.64305, 4.30333, 2.42879] * u.km / u.s
    ss0 = State.from_vectors(Earth, r0, v0)
    tof = 0 * u.s
    or0 = Orbit(ss0)

    or1 = or0.propagate(tof)

    r, v = or1.state.rv()

    assert_array_almost_equal(r.value, r0.value)
    assert_array_almost_equal(v.value, v0.value)


def test_apply_zero_maneuver_returns_equal_state():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    ss = State.from_classical(Sun, _d, _, _a, _a, _a, _a)
    dt = 0 * u.s
    dv = [0, 0, 0] * u.km / u.s
    orbit = Orbit(ss)
    orbit_new = orbit.apply_maneuver([(dt, dv)])
    assert_almost_equal(orbit_new.state.r.to(u.km).value,
                        ss.r.to(u.km).value)
    assert_almost_equal(orbit_new.state.v.to(u.km / u.s).value,
                        ss.v.to(u.km / u.s).value)
