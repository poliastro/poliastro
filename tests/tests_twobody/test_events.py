import numpy as np
import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time
from numpy.linalg import norm

from poliastro.bodies import Earth
from poliastro.constants import H0_earth, rho0_earth
from poliastro.core.perturbations import atmospheric_drag_exponential
from poliastro.core.propagation import func_twobody
from poliastro.twobody import Orbit
from poliastro.twobody.events import (
    AltitudeCrossEvent,
    LatitudeCrossEvent,
    NodeCrossEvent,
    PenumbraEvent,
    UmbraEvent,
)
from poliastro.twobody.propagation import cowell


@pytest.mark.slow
def test_altitude_crossing():
    # Test decreasing altitude cross over Earth. No analytic solution.
    R = Earth.R.to(u.km).value

    orbit = Orbit.circular(Earth, 230 * u.km)
    t_flight = 48.209538 * u.d

    # Parameters of a body
    C_D = 2.2  # Dimensionless (any value would do)
    A_over_m = ((np.pi / 4.0) * (u.m ** 2) / (100 * u.kg)).to_value(
        u.km ** 2 / u.kg
    )  # km^2/kg

    # Parameters of the atmosphere
    rho0 = rho0_earth.to(u.kg / u.km ** 3).value  # kg/km^3
    H0 = H0_earth.to(u.km).value  # km

    tofs = [50] * u.d

    thresh_alt = 50  # km
    altitude_cross_event = AltitudeCrossEvent(thresh_alt, R)
    events = [altitude_cross_event]

    def f(t0, u_, k):
        du_kep = func_twobody(t0, u_, k)
        ax, ay, az = atmospheric_drag_exponential(
            t0, u_, k, R=R, C_D=C_D, A_over_m=A_over_m, H0=H0, rho0=rho0
        )
        du_ad = np.array([0, 0, 0, ax, ay, az])
        return du_kep + du_ad

    rr, _ = cowell(
        Earth.k,
        orbit.r,
        orbit.v,
        tofs,
        events=events,
        f=f,
    )

    assert_quantity_allclose(norm(rr[0].to(u.km).value) - thresh_alt, R)
    assert_quantity_allclose(altitude_cross_event.last_t, t_flight, rtol=1e-2)


def test_altitude_cross_not_happening_is_ok():
    R = Earth.R.to(u.km).value

    orbit = Orbit.circular(Earth, 230 * u.km)

    tofs = [25] * u.d

    thresh_alt = 50  # km
    altitude_cross_event = AltitudeCrossEvent(thresh_alt, R)
    events = [altitude_cross_event]

    rr, _ = cowell(
        Earth.k,
        orbit.r,
        orbit.v,
        tofs,
        events=events,
    )

    assert altitude_cross_event.last_t == tofs[-1]


def test_latitude_cross_event():
    r = [-6142438.668, 3492467.56, -25767.257] << u.km
    v = [505.848, 942.781, 7435.922] << u.km / u.s
    orbit = Orbit.from_vectors(Earth, r, v)

    thresh_lat = 60 * u.deg
    latitude_cross_event = LatitudeCrossEvent(orbit, thresh_lat, terminal=True)
    t_lat = 1701.716842130476 * u.s

    tofs = [5] * u.d

    events = [latitude_cross_event]
    rr, _ = cowell(
        Earth.k,
        orbit.r,
        orbit.v,
        tofs,
        events=events,
    )

    assert_quantity_allclose(latitude_cross_event.last_t, t_lat)


def test_penumbra_event_not_triggering_is_ok():
    attractor = Earth
    tof = 100 * u.s
    r0 = np.array([281.89, 1411.473, 750.672])
    v0 = np.array([7.36138, 2.98997, 1.64354])
    orbit = Orbit.from_vectors(attractor, r0 * u.km, v0 * u.km / u.s)

    penumbra_event = PenumbraEvent(orbit)
    events = [penumbra_event]

    rr, _ = cowell(
        attractor.k,
        orbit.r,
        orbit.v,
        [tof] * u.s,
        events=events,
    )

    assert penumbra_event.last_t == tof


def test_umbra_event_not_triggering_is_ok():
    attractor = Earth
    tof = 100 * u.s
    r0 = np.array([281.89, 1411.473, 750.672])
    v0 = np.array([7.36138, 2.98997, 1.64354])
    orbit = Orbit.from_vectors(attractor, r0 * u.km, v0 * u.km / u.s)

    umbra_event = UmbraEvent(orbit)
    events = [umbra_event]

    rr, _ = cowell(
        attractor.k,
        orbit.r,
        orbit.v,
        [tof] * u.s,
        events=events,
    )

    assert umbra_event.last_t == tof


def test_umbra_event_crossing():
    expected_umbra_t = Time("2020-01-01 00:04:51.328", scale="utc")  # From Orekit.
    attractor = Earth
    tof = 2 * u.d
    epoch = Time("2020-01-01", scale="utc")
    coe = (
        6828137.0 * u.m,
        0.0073 * u.one,
        87.0 * u.deg,
        20.0 * u.deg,
        10.0 * u.deg,
        0 * u.deg,
    )
    orbit = Orbit.from_classical(attractor, *coe, epoch=epoch)

    umbra_event = UmbraEvent(orbit, terminal=True)
    events = [umbra_event]

    rr, _ = cowell(
        attractor.k,
        orbit.r,
        orbit.v,
        [tof] * u.s,
        events=events,
    )

    assert expected_umbra_t.isclose(epoch + umbra_event.last_t, atol=1 * u.s)


def test_penumbra_event_crossing():
    expected_penumbra_t = Time("2020-01-01 00:04:26.060", scale="utc")  # From Orekit.
    attractor = Earth
    tof = 2 * u.d
    epoch = Time("2020-01-01", scale="utc")
    coe = (
        6828137.0 * u.m,
        0.0073 * u.one,
        87.0 * u.deg,
        20.0 * u.deg,
        10.0 * u.deg,
        0 * u.deg,
    )
    orbit = Orbit.from_classical(attractor, *coe, epoch=epoch)

    penumbra_event = PenumbraEvent(orbit, terminal=True)
    events = [penumbra_event]

    rr, _ = cowell(
        attractor.k,
        orbit.r,
        orbit.v,
        [tof] * u.s,
        events=events,
    )

    assert expected_penumbra_t.isclose(epoch + penumbra_event.last_t, atol=1 * u.s)


def test_node_cross_event():
    t_node = 3.46524036 * u.s
    r = [-6142438.668, 3492467.56, -25767.257] << u.km
    v = [505.848, 942.781, 7435.922] << u.km / u.s
    orbit = Orbit.from_vectors(Earth, r, v)

    node_event = NodeCrossEvent(terminal=True)
    events = [node_event]

    tofs = [0.01, 0.1, 0.5, 0.8, 1, 3, 5, 6, 10, 15] << u.s
    rr, vv = cowell(
        Earth.k,
        orbit.r,
        orbit.v,
        tofs,
        events=events,
    )

    assert_quantity_allclose(node_event.last_t, t_node)


def test_node_event_equatorial_orbit():
    node_event = NodeCrossEvent(terminal=True)
    events = [node_event]

    r = np.array([9946.2, 1035.4, 0.0])
    v = np.array([7.0, -0.1, 0.0])
    orb = Orbit.from_vectors(Earth, r * u.km, v * u.km / u.s)

    tofs = [5, 10, 50] << u.s
    rr, vv = cowell(
        Earth.k,
        orb.r,
        orb.v,
        tofs,
        events=events,
    )

    assert_quantity_allclose(node_event.last_t, 0.0 * u.s, atol=1e-1 * u.s)
