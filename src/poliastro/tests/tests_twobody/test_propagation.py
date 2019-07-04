import numpy as np
import pytest
from astropy import time, units as u
from astropy.coordinates import CartesianRepresentation
from astropy.tests.helper import assert_quantity_allclose
from numpy.testing import assert_allclose
from pytest import approx

from poliastro.bodies import Earth, Moon, Sun
from poliastro.constants import J2000
from poliastro.core.elements import rv2coe
from poliastro.examples import iss
from poliastro.frames import Planes
from poliastro.twobody import Orbit
from poliastro.twobody.propagation import (
    cowell,
    danby,
    goodingI,
    kepler,
    markley,
    mean_motion,
    mikkola,
    pimienta,
)
from poliastro.util import norm

propagators_set = [mean_motion, kepler, mikkola, markley, pimienta, goodingI, danby]
accurate_propagators_elliptic = [markley, pimienta, goodingI, danby]
accurate_propagators_hyperbolic = [kepler, danby]


@pytest.mark.parametrize("ecc", [0.99, 0.99, 0.999, 0.9999])
@pytest.mark.parametrize("propagator", accurate_propagators_elliptic)
def test_elliptic_near_parabolic(ecc, propagator):
    _a = 0.0 * u.rad
    tof = 1.0 * u.min
    ss0 = Orbit.from_classical(
        Earth, 10000 * u.km, ecc * u.one, _a, _a, _a, 1.0 * u.rad
    )

    ss_cowell = ss0.propagate(tof, method=cowell)
    ss_propagator = ss0.propagate(tof, method=propagator)

    assert_quantity_allclose(ss_propagator.r, ss_cowell.r)
    assert_quantity_allclose(ss_propagator.v, ss_cowell.v)


@pytest.mark.parametrize("ecc", [1.0001, 1.001, 1.01, 1.1])
@pytest.mark.parametrize("propagator", accurate_propagators_hyperbolic)
def test_hyperbolic_near_parabolic(ecc, propagator):
    _a = 0.0 * u.rad
    tof = 1.0 * u.min
    ss0 = Orbit.from_classical(
        Earth, -10000 * u.km, ecc * u.one, _a, _a, _a, 1.0 * u.rad
    )

    ss_cowell = ss0.propagate(tof, method=cowell)
    ss_propagator = ss0.propagate(tof, method=propagator)

    assert_quantity_allclose(ss_propagator.r, ss_cowell.r)
    assert_quantity_allclose(ss_propagator.v, ss_cowell.v)


@pytest.mark.parametrize("propagator", [markley, pimienta])
def test_near_equatorial(propagator):
    r = [8.0e3, 1.0e3, 0.0] * u.km
    v = [-0.5, -0.5, 0.0001] * u.km / u.s
    tof = 1.0 * u.h
    ss0 = Orbit.from_vectors(Earth, r, v)

    ss_cowell = ss0.propagate(tof, method=cowell)
    ss_propagator = ss0.propagate(tof, method=propagator)

    assert_quantity_allclose(ss_propagator.r, ss_cowell.r, rtol=1e-4)
    assert_quantity_allclose(ss_propagator.v, ss_cowell.v, rtol=1e-4)


@pytest.mark.parametrize("propagator", propagators_set)
def test_propagation(propagator):
    # Data from Vallado, example 2.4
    r0 = [1131.340, -2282.343, 6672.423] * u.km
    v0 = [-5.64305, 4.30333, 2.42879] * u.km / u.s
    expected_r = [-4219.7527, 4363.0292, -3958.7666] * u.km
    expected_v = [3.689866, -1.916735, -6.112511] * u.km / u.s

    ss0 = Orbit.from_vectors(Earth, r0, v0)
    tof = 40 * u.min
    ss1 = ss0.propagate(tof, method=propagator)

    r, v = ss1.rv()

    assert_quantity_allclose(r, expected_r, rtol=1e-5)
    assert_quantity_allclose(v, expected_v, rtol=1e-4)


def test_propagating_to_certain_nu_is_correct():
    # take an elliptic orbit
    a = 1.0 * u.AU
    ecc = 1.0 / 3.0 * u.one
    _a = 0.0 * u.rad
    nu = 10 * u.deg
    elliptic = Orbit.from_classical(Sun, a, ecc, _a, _a, _a, nu)

    elliptic_at_perihelion = elliptic.propagate_to_anomaly(0.0 * u.rad)
    r_per, _ = elliptic_at_perihelion.rv()

    elliptic_at_aphelion = elliptic.propagate_to_anomaly(np.pi * u.rad)
    r_ap, _ = elliptic_at_aphelion.rv()

    assert_quantity_allclose(norm(r_per), a * (1.0 - ecc))
    assert_quantity_allclose(norm(r_ap), a * (1.0 + ecc))

    # TODO: Test specific values
    assert elliptic_at_perihelion.epoch > elliptic.epoch
    assert elliptic_at_aphelion.epoch > elliptic.epoch

    # test 10 random true anomaly values
    for _ in range(10):
        nu = np.random.uniform(low=0.0, high=2 * np.pi)
        elliptic = elliptic.propagate_to_anomaly(nu * u.rad)
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


@pytest.mark.xfail
def test_propagation_mean_motion_parabolic_fails():
    # example from Howard Curtis (3rd edition), section 3.5, problem 3.15
    p = 2.0 * 6600 * u.km
    _a = 0.0 * u.deg
    orbit = Orbit.parabolic(Earth, p, _a, _a, _a, _a)
    orbit = orbit.propagate(0.8897 / 2.0 * u.h, method=mean_motion)

    _, _, _, _, _, nu0 = rv2coe(
        Earth.k.to(u.km ** 3 / u.s ** 2).value,
        orbit.r.to(u.km).value,
        orbit.v.to(u.km / u.s).value,
    )
    assert_quantity_allclose(nu0, np.deg2rad(90.0), rtol=1e-4)

    orbit = Orbit.parabolic(Earth, p, _a, _a, _a, _a)
    orbit = orbit.propagate(36.0 * u.h, method=mean_motion)
    assert_quantity_allclose(norm(orbit.r), 304700.0 * u.km, rtol=1e-4)


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


def test_propagation_hyperbolic_zero_time_returns_same_state():
    ss0 = Orbit.from_classical(
        Earth,
        -27112.5464 * u.km,
        1.25 * u.one,
        0 * u.deg,
        0 * u.deg,
        0 * u.deg,
        0 * u.deg,
    )
    r0, v0 = ss0.rv()
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
    assert_allclose(orbit_new.r.to(u.km).value, ss.r.to(u.km).value)
    assert_allclose(orbit_new.v.to(u.km / u.s).value, ss.v.to(u.km / u.s).value)


def test_cowell_propagation_with_zero_acceleration_equals_kepler():
    # Data from Vallado, example 2.4

    r0 = np.array([1131.340, -2282.343, 6672.423]) * u.km
    v0 = np.array([-5.64305, 4.30333, 2.42879]) * u.km / u.s
    tofs = [40 * 60.0] * u.s

    orbit = Orbit.from_vectors(Earth, r0, v0)

    expected_r = np.array([-4219.7527, 4363.0292, -3958.7666]) * u.km
    expected_v = np.array([3.689866, -1.916735, -6.112511]) * u.km / u.s

    r, v = cowell(Earth.k, orbit.r, orbit.v, tofs, ad=None)

    assert_quantity_allclose(r[0], expected_r, rtol=1e-5)
    assert_quantity_allclose(v[0], expected_v, rtol=1e-4)


def test_cowell_propagation_circle_to_circle():
    # From [Edelbaum, 1961]
    accel = 1e-7

    def constant_accel(t0, u, k):
        v = u[3:]
        norm_v = (v[0] ** 2 + v[1] ** 2 + v[2] ** 2) ** 0.5
        return accel * v / norm_v

    ss = Orbit.circular(Earth, 500 * u.km)
    tofs = [20] * ss.period

    r, v = cowell(Earth.k, ss.r, ss.v, tofs, ad=constant_accel)

    ss_final = Orbit.from_vectors(Earth, r[0], v[0])

    da_a0 = (ss_final.a - ss.a) / ss.a
    dv_v0 = abs(norm(ss_final.v) - norm(ss.v)) / norm(ss.v)
    assert_quantity_allclose(da_a0, 2 * dv_v0, rtol=1e-2)

    dv = abs(norm(ss_final.v) - norm(ss.v))
    accel_dt = accel * u.km / u.s ** 2 * tofs[0]
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


# @pytest.mark.filterwarnings("ignore:ERFA")
# @pytest.mark.filterwarnings("ignore::UserWarning")
# @pytest.mark.parametrize(
#    "method",
#    [
#        mean_motion,
#        pytest.param(kepler, marks=pytest.mark.xfail),
#        pytest.param(
#            cowell,
#            marks=pytest.mark.skip(reason="for now propagation takes too long time"),
#        ),
#    ],
# )
@pytest.mark.parametrize("propagator", [danby, markley, goodingI])
def test_propagate_long_times_keeps_geometry(propagator):
    # See https://github.com/poliastro/poliastro/issues/265
    time_of_flight = 100 * u.year

    res = iss.propagate(time_of_flight, method=propagator)

    assert_quantity_allclose(iss.a, res.a)
    assert_quantity_allclose(iss.ecc, res.ecc)
    assert_quantity_allclose(iss.inc, res.inc)
    assert_quantity_allclose(iss.raan, res.raan)
    assert_quantity_allclose(iss.argp, res.argp)

    assert_quantity_allclose(
        (res.epoch - iss.epoch).to(time_of_flight.unit), time_of_flight
    )


@pytest.mark.filterwarnings("ignore::astropy._erfa.core.ErfaWarning")
def test_long_propagations_kepler_agrees_mean_motion():
    tof = 100 * u.year
    r_mm, v_mm = iss.propagate(tof, method=mean_motion).rv()
    r_k, v_k = iss.propagate(tof, method=kepler).rv()
    assert_quantity_allclose(r_mm, r_k)
    assert_quantity_allclose(v_mm, v_k)

    r_halleys = [-9018878.63569932, -94116054.79839276, 22619058.69943215]  # km
    v_halleys = [-49.95092305, -12.94843055, -4.29251577]  # km/s
    halleys = Orbit.from_vectors(Sun, r_halleys * u.km, v_halleys * u.km / u.s)

    r_mm, v_mm = halleys.propagate(tof, method=mean_motion).rv()
    r_k, v_k = halleys.propagate(tof, method=kepler).rv()
    assert_quantity_allclose(r_mm, r_k)
    assert_quantity_allclose(v_mm, v_k)


@pytest.mark.parametrize("method", [mean_motion, kepler])
def test_long_propagation_preserves_orbit_elements(method):
    tof = 100 * u.year
    r_halleys = np.array(
        [-9018878.63569932, -94116054.79839276, 22619058.69943215]
    )  # km
    v_halleys = np.array([-49.95092305, -12.94843055, -4.29251577])  # km/s
    halleys = Orbit.from_vectors(Sun, r_halleys * u.km, v_halleys * u.km / u.s)

    params_ini = rv2coe(Sun.k.to(u.km ** 3 / u.s ** 2).value, r_halleys, v_halleys)[:-1]
    r_new, v_new = halleys.propagate(tof, method=method).rv()
    params_final = rv2coe(
        Sun.k.to(u.km ** 3 / u.s ** 2).value,
        r_new.to(u.km).value,
        v_new.to(u.km / u.s).value,
    )[:-1]
    assert_quantity_allclose(params_ini, params_final)


def test_propagation_sets_proper_epoch():
    expected_epoch = time.Time("2017-09-01 12:05:50", scale="tdb")

    r = [-2.76132873e08, -1.71570015e08, -1.09377634e08] * u.km
    v = [13.17478674, -9.82584125, -1.48126639] * u.km / u.s
    florence = Orbit.from_vectors(Sun, r, v, plane=Planes.EARTH_ECLIPTIC)

    propagated = florence.propagate(expected_epoch)

    assert propagated.epoch == expected_epoch


def test_sample_custom_body_raises_warning_and_returns_coords():
    # See https://github.com/poliastro/poliastro/issues/649
    orbit = Orbit.circular(Moon, 100 * u.km)

    with pytest.warns(
        UserWarning, match="returning only cartesian coordinates instead"
    ):
        coords = orbit.sample(10)

    assert isinstance(coords, CartesianRepresentation)
    assert len(coords) == 10


def test_propagation_custom_body_works():
    # See https://github.com/poliastro/poliastro/issues/649
    orbit = Orbit.circular(Moon, 100 * u.km)
    orbit.propagate(1 * u.h)
