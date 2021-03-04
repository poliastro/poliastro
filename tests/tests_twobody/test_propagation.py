import numpy as np
import pytest
from astropy import time, units as u
from astropy.coordinates import CartesianRepresentation
from astropy.tests.helper import assert_quantity_allclose
from hypothesis import given, settings, strategies as st
from numpy.testing import assert_allclose
from pytest import approx

from poliastro.bodies import Earth, Moon, Sun
from poliastro.constants import J2000
from poliastro.core.elements import rv2coe
from poliastro.core.propagation import (
    danby_coe,
    func_twobody,
    gooding_coe,
    markley_coe,
    mikkola_coe,
    pimienta_coe,
)
from poliastro.core.propagation.farnocchia import farnocchia_coe
from poliastro.examples import iss
from poliastro.frames import Planes
from poliastro.twobody import Orbit
from poliastro.twobody.propagation import (
    ALL_PROPAGATORS,
    ELLIPTIC_PROPAGATORS,
    HYPERBOLIC_PROPAGATORS,
    PARABOLIC_PROPAGATORS,
    cowell,
    danby,
    danby_classical,
    farnocchia,
    farnocchia_classical,
    gooding,
    gooding_classical,
    markley,
    markley_classical,
    mikkola,
    mikkola_classical,
    pimienta,
    pimienta_classical,
    vallado,
)
from poliastro.util import norm


@pytest.fixture(scope="module")
def halley():
    return Orbit.from_vectors(
        Sun,
        [-9018878.63569932, -94116054.79839276, 22619058.69943215] * u.km,
        [-49.95092305, -12.94843055, -4.29251577] * u.km / u.s,
    )


@pytest.mark.parametrize("ecc", [0.9, 0.99, 0.999, 0.9999, 0.99999])
@pytest.mark.parametrize("propagator", ELLIPTIC_PROPAGATORS)
def test_elliptic_near_parabolic(ecc, propagator):
    # 'kepler fails if really close to parabolic'. Refer to issue #714.
    if propagator in [vallado] and ecc > 0.99:
        pytest.xfail()

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
@pytest.mark.parametrize("propagator", HYPERBOLIC_PROPAGATORS)
def test_hyperbolic_near_parabolic(ecc, propagator):
    # Still not implemented. Refer to issue #714.
    if propagator in [pimienta, gooding, gooding_classical]:
        pytest.skip()

    _a = 0.0 * u.rad
    tof = 1.0 * u.min
    ss0 = Orbit.from_classical(
        Earth, -10000 * u.km, ecc * u.one, _a, _a, _a, 1.0 * u.rad
    )

    ss_cowell = ss0.propagate(tof, method=cowell)
    ss_propagator = ss0.propagate(tof, method=propagator)

    assert_quantity_allclose(ss_propagator.r, ss_cowell.r)
    assert_quantity_allclose(ss_propagator.v, ss_cowell.v)


@pytest.mark.parametrize("propagator", [markley])
def test_near_equatorial(propagator):
    r = [8.0e3, 1.0e3, 0.0] * u.km
    v = [-0.5, -0.5, 0.0001] * u.km / u.s
    tof = 1.0 * u.h
    ss0 = Orbit.from_vectors(Earth, r, v)

    ss_cowell = ss0.propagate(tof, method=cowell)
    ss_propagator = ss0.propagate(tof, method=propagator)

    assert_quantity_allclose(ss_propagator.r, ss_cowell.r, rtol=1e-4)
    assert_quantity_allclose(ss_propagator.v, ss_cowell.v, rtol=1e-4)


@pytest.mark.parametrize("propagator", ALL_PROPAGATORS)
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


def test_propagation_cartesian_and_classical_matches():
    propagators = [farnocchia, mikkola, markley, pimienta, gooding, danby]
    propagators_classical = [
        farnocchia_classical,
        mikkola_classical,
        markley_classical,
        pimienta_classical,
        gooding_classical,
        danby_classical,
    ]
    # Data from Vallado, example 2.4
    r0 = [1131.340, -2282.343, 6672.423] * u.km
    v0 = [-5.64305, 4.30333, 2.42879] * u.km / u.s

    ss0 = Orbit.from_vectors(Earth, r0, v0)
    tof = 40 * u.min

    for propagator_pair in zip(propagators, propagators_classical):
        ss1 = ss0.propagate(tof, method=propagator_pair[0])
        ss2 = ss0.propagate(tof, method=propagator_pair[1])

        r1, v1 = ss1.rv()
        r2, v2 = ss2.rv()

        assert_quantity_allclose(r1, r2, rtol=1e-2)
        assert_quantity_allclose(v1, v2, rtol=1e-2)


def test_propagating_to_certain_nu_is_correct():
    # Take an elliptic orbit
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

    # Test 10 random true anomaly values
    # TODO: Rework this test
    for nu in np.random.uniform(low=-np.pi, high=np.pi, size=10):
        elliptic = elliptic.propagate_to_anomaly(nu * u.rad)
        r, _ = elliptic.rv()
        assert_quantity_allclose(norm(r), a * (1.0 - ecc ** 2) / (1 + ecc * np.cos(nu)))


def test_propagate_to_anomaly_in_the_past_fails_for_open_orbits():
    r0 = [Earth.R.to(u.km).value + 300, 0, 0] * u.km
    v0 = [0, 15, 0] * u.km / u.s
    orb = Orbit.from_vectors(Earth, r0, v0)

    with pytest.raises(ValueError, match="True anomaly -0.02 rad not reachable"):
        orb.propagate_to_anomaly(orb.nu - 1 * u.deg)


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


@pytest.mark.parametrize("propagator", PARABOLIC_PROPAGATORS)
def test_propagation_parabolic(propagator):
    # Example from Howard Curtis (3rd edition), section 3.5, problem 3.15
    # TODO: add parabolic solver in some parabolic propagators, refer to #417
    if propagator in [mikkola, gooding, gooding_classical]:
        pytest.skip()

    p = 2.0 * 6600 * u.km
    _a = 0.0 * u.deg
    orbit = Orbit.parabolic(Earth, p, _a, _a, _a, _a)
    orbit = orbit.propagate(0.8897 / 2.0 * u.h, method=propagator)

    _, _, _, _, _, nu0 = rv2coe(
        Earth.k.to(u.km ** 3 / u.s ** 2).value,
        orbit.r.to(u.km).value,
        orbit.v.to(u.km / u.s).value,
    )
    assert_quantity_allclose(nu0, np.deg2rad(90.0), rtol=1e-4)

    orbit = Orbit.parabolic(Earth, p, _a, _a, _a, _a)
    orbit = orbit.propagate(36.0 * u.h, method=propagator)
    assert_quantity_allclose(norm(orbit.r), 304700.0 * u.km, rtol=1e-4)


def test_propagation_zero_time_returns_same_state():
    # Bug #50
    r0 = [1131.340, -2282.343, 6672.423] * u.km  # type: u.Quantity
    v0 = [-5.64305, 4.30333, 2.42879] * u.km / u.s
    ss0 = Orbit.from_vectors(Earth, r0, v0)
    tof = 0 * u.s

    ss1 = ss0.propagate(tof)

    r, v = ss1.rv()

    assert_quantity_allclose(r, r0)
    assert_quantity_allclose(v, v0)


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

    assert_quantity_allclose(r, r0)
    assert_quantity_allclose(v, v0)


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

    r, v = cowell(Earth.k, orbit.r, orbit.v, tofs)

    assert_quantity_allclose(r[0], expected_r, rtol=1e-5)
    assert_quantity_allclose(v[0], expected_v, rtol=1e-4)


def test_cowell_propagation_circle_to_circle():
    # From [Edelbaum, 1961]
    accel = 1e-7

    def constant_accel(t0, u_, k):
        v = u_[3:]
        norm_v = (v[0] ** 2 + v[1] ** 2 + v[2] ** 2) ** 0.5
        return accel * v / norm_v

    def f(t0, u_, k):
        du_kep = func_twobody(t0, u_, k)
        ax, ay, az = constant_accel(t0, u_, k)
        du_ad = np.array([0, 0, 0, ax, ay, az])

        return du_kep + du_ad

    ss = Orbit.circular(Earth, 500 * u.km)
    tofs = [20] * ss.period

    r, v = cowell(Earth.k, ss.r, ss.v, tofs, f=f)

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


@pytest.mark.parametrize("propagator", [danby, markley, gooding])
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


@pytest.mark.filterwarnings("ignore::erfa.core.ErfaWarning")
def test_long_propagations_vallado_agrees_farnocchia():
    tof = 100 * u.year
    r_mm, v_mm = iss.propagate(tof, method=farnocchia).rv()
    r_k, v_k = iss.propagate(tof, method=vallado).rv()
    assert_quantity_allclose(r_mm, r_k)
    assert_quantity_allclose(v_mm, v_k)

    r_halleys = [-9018878.63569932, -94116054.79839276, 22619058.69943215]  # km
    v_halleys = [-49.95092305, -12.94843055, -4.29251577]  # km/s
    halleys = Orbit.from_vectors(Sun, r_halleys * u.km, v_halleys * u.km / u.s)

    r_mm, v_mm = halleys.propagate(tof, method=farnocchia).rv()
    r_k, v_k = halleys.propagate(tof, method=vallado).rv()
    assert_quantity_allclose(r_mm, r_k)
    assert_quantity_allclose(v_mm, v_k)


@st.composite
def with_units(draw, elements, unit):
    value = draw(elements)
    return value * unit


@settings(deadline=None)
@given(
    tof=with_units(
        elements=st.floats(
            min_value=80, max_value=120, allow_nan=False, allow_infinity=False
        ),
        unit=u.year,
    )
)
@pytest.mark.parametrize("method", [farnocchia, vallado])
def test_long_propagation_preserves_orbit_elements(tof, method, halley):
    expected_slow_classical = halley.classical()[:-1]

    slow_classical = halley.propagate(tof, method=method).classical()[:-1]

    for element, expected_element in zip(slow_classical, expected_slow_classical):
        assert_quantity_allclose(element, expected_element)


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

    coords = orbit.sample(10)

    assert isinstance(coords, CartesianRepresentation)
    assert len(coords) == 10


def test_propagation_custom_body_works():
    # See https://github.com/poliastro/poliastro/issues/649
    orbit = Orbit.circular(Moon, 100 * u.km)
    orbit.propagate(1 * u.h)


@pytest.mark.parametrize(
    "propagator_coe",
    [danby_coe, markley_coe, pimienta_coe, mikkola_coe, farnocchia_coe, gooding_coe],
)
def test_propagate_with_coe(propagator_coe):

    period = iss.period
    a, ecc, inc, raan, argp, nu = iss.classical()
    p = a * (1 - ecc ** 2)

    # Delete the units
    p = p.to_value(u.km)
    ecc = ecc.value
    period = period.to_value(u.s)
    inc = inc.to_value(u.rad)
    raan = raan.to_value(u.rad)
    argp = argp.to_value(u.rad)
    nu = nu.to_value(u.rad)
    k = iss.attractor.k.to_value(u.km ** 3 / u.s ** 2)

    r0, v0 = propagator_coe(k, p, ecc, inc, raan, argp, nu, period)
    _, _, _, _, _, nu_final = rv2coe(k, r0, v0)

    assert_quantity_allclose(nu_final, nu)
