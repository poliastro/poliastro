import numpy as np
import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time
from numpy.linalg import norm

from poliastro.bodies import Earth
from poliastro.constants import H0_earth, rho0_earth
from poliastro.core.events import line_of_sight
from poliastro.core.perturbations import atmospheric_drag_exponential
from poliastro.core.propagation import func_twobody
from poliastro.twobody import Orbit
from poliastro.twobody.events import (
    AltitudeCrossEvent,
    LatitudeCrossEvent,
    LithobrakeEvent,
    LosEvent,
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
    orbit = Orbit.from_classical(
        attractor=attractor,
        a=6828137.0 * u.m,
        ecc=0.0073 * u.one,
        inc=87.0 * u.deg,
        raan=20.0 * u.deg,
        argp=10.0 * u.deg,
        nu=0 * u.deg,
        epoch=epoch
    )

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
    orbit = Orbit.from_classical(
        attractor=attractor,
        a=6828137.0 * u.m,
        ecc=0.0073 * u.one,
        inc=87.0 * u.deg,
        raan=20.0 * u.deg,
        argp=10.0 * u.deg,
        nu=0 * u.deg,
        epoch=epoch
    )

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


def test_orbit_propagation_continues_if_events_terminal_is_False():
    r = [-6142438.668, 3492467.56, -25767.257] << u.km
    v = [505.848, 942.781, 7435.922] << u.km / u.s
    orbit = Orbit.from_vectors(Earth, r, v)

    thresh_lat = 60 * u.deg
    # Event occurs at ~1701.7 s.
    latitude_cross_event = LatitudeCrossEvent(orbit, thresh_lat, terminal=False)
    events = [latitude_cross_event]

    # The last two tofs are after the detection of the event.
    tofs = [1000, 1250, 1500, 1710, 2000] << u.s
    rr, _ = cowell(
        Earth.k,
        orbit.r,
        orbit.v,
        tofs,
        events=events,
    )

    # Check position vectors don't repeat during propagation.
    assert not np.allclose(rr[-1], rr[-2])


def test_orbit_propagation_position_vector_does_not_repeat_if_events_terminal_is_True():
    r = [-6142438.668, 3492467.56, -25767.257] << u.km
    v = [505.848, 942.781, 7435.922] << u.km / u.s
    orbit = Orbit.from_vectors(Earth, r, v)

    thresh_lat = 60 * u.deg
    # Event occurs at ~1701.7 s.
    latitude_cross_event = LatitudeCrossEvent(orbit, thresh_lat, terminal=True)
    events = [latitude_cross_event]

    # The last two tofs are after the detection of the event.
    tofs = [1000, 1250, 1500, 1710, 2000] << u.s
    rr, _ = cowell(
        Earth.k,
        orbit.r,
        orbit.v,
        tofs,
        events=events,
    )

    # Check position vector doesn't repeat if terminal set to True.
    assert len(rr) == 4  # From the 5th tof in tofs, position vector starts repeating.


@pytest.mark.parametrize(
    "latitude_terminal,penumbra_terminal,rr_length,t_end",
    [
        (True, True, 4, 266.15058 * u.s),
        (True, False, 5, 305.65173 * u.s),
        (False, True, 4, 266.15058 * u.s),
        (False, False, 6, 500 * u.s),
    ],
)
def test_propagation_stops_if_atleast_one_event_has_terminal_set_to_True(
    latitude_terminal, penumbra_terminal, rr_length, t_end
):
    # Penumbra occurs at 266.15058s and latitude event occurs at 305.65173s.
    # `terminals` is for latitude event and penumbra event, in that order.
    attractor = Earth
    tofs = [50, 100, 150, 300, 400, 500] << u.s
    epoch = Time("2020-01-01", scale="utc")
    orbit = Orbit.from_classical(
        attractor=attractor,
        a=6828137.0 * u.m,
        ecc=0.0073 * u.one,
        inc=87.0 * u.deg,
        raan=20.0 * u.deg,
        argp=10.0 * u.deg,
        nu=0 * u.deg,
        epoch=epoch
    )

    penumbra_event = PenumbraEvent(orbit, terminal=penumbra_terminal)

    thresh_lat = 30 * u.deg
    latitude_cross_event = LatitudeCrossEvent(
        orbit, thresh_lat, terminal=latitude_terminal
    )
    events = [penumbra_event, latitude_cross_event]

    rr, _ = cowell(
        attractor.k,
        orbit.r,
        orbit.v,
        tofs,
        events=events,
    )

    assert len(rr) == rr_length
    if penumbra_terminal:
        assert_quantity_allclose(penumbra_event.last_t, t_end)
    elif latitude_terminal and not penumbra_terminal:
        assert_quantity_allclose(latitude_cross_event.last_t, t_end)
    else:
        assert_quantity_allclose(t_end, tofs[-1])


def test_line_of_sight():
    # From Vallado example 5.6
    r1 = np.array([0, -4464.696, -5102.509]) << u.km
    r2 = np.array([0, 5740.323, 3189.068]) << u.km
    r_sun = np.array([122233179, -76150708, 33016374]) << u.km
    R = Earth.R.to(u.km).value

    los = line_of_sight(r1.value, r2.value, R)
    los_with_sun = line_of_sight(r1.value, r_sun.value, R)

    assert los < 0  # No LOS condition.
    assert los_with_sun >= 0  # LOS condition.


def test_LOS_event_raises_warning_if_norm_of_r1_less_than_attractor_radius_during_propagation():
    r2 = np.array([-500, 1500, 4012.09]) << u.km
    v2 = np.array([5021.38, -2900.7, 1000.354]) << u.km / u.s
    orbit = Orbit.from_vectors(Earth, r2, v2)

    tofs = [100, 500, 1000, 2000] << u.s
    # Propagate the secondary body to generate its position coordinates.
    rr, vv = cowell(
        Earth.k,
        orbit.r,
        orbit.v,
        tofs,
    )
    pos_coords = rr  # Trajectory of the secondary body.

    r1 = (
        np.array([0, -5010.696, -5102.509]) << u.km
    )  # This position vectors' norm gets less than attractor radius.
    v1 = np.array([736.138, 29899.7, 164.354]) << u.km / u.s
    orb = Orbit.from_vectors(Earth, r1, v1)

    los_event = LosEvent(Earth, pos_coords, terminal=True)
    events = [los_event]
    tofs = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.5] << u.s

    with pytest.warns(UserWarning, match="The norm of the position vector"):
        r, v = cowell(
            Earth.k,
            orb.r,
            orb.v,
            tofs,
            events=events,
        )


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_LOS_event_with_lithobrake_event_raises_warning_when_satellite_cuts_attractor():
    r2 = np.array([-500, 1500, 4012.09]) << u.km
    v2 = np.array([5021.38, -2900.7, 1000.354]) << u.km / u.s
    orbit = Orbit.from_vectors(Earth, r2, v2)

    tofs = [100, 500, 1000, 2000] << u.s
    # Propagate the secondary body to generate its position coordinates.
    rr, vv = cowell(
        Earth.k,
        orbit.r,
        orbit.v,
        tofs,
    )
    pos_coords = rr  # Trajectory of the secondary body.

    r1 = np.array([0, -5010.696, -5102.509]) << u.km
    v1 = np.array([736.138, 2989.7, 164.354]) << u.km / u.s
    orb = Orbit.from_vectors(Earth, r1, v1)

    los_event = LosEvent(Earth, pos_coords, terminal=True)
    tofs = [0.003, 0.004, 0.01, 0.02, 0.03, 0.04, 0.07, 0.1, 0.2, 0.3, 0.4, 1, 3] << u.s

    lithobrake_event = LithobrakeEvent(Earth.R.to_value(u.km))
    events = [lithobrake_event, los_event]
    r, v = cowell(
        Earth.k,
        orb.r,
        orb.v,
        tofs,
        events=events,
    )

    assert lithobrake_event.last_t < los_event.last_t


def test_LOS_event():
    t_los = 2327.165 * u.s
    r2 = np.array([-500, 1500, 4012.09]) << u.km
    v2 = np.array([5021.38, -2900.7, 1000.354]) << u.km / u.s
    orbit = Orbit.from_vectors(Earth, r2, v2)

    tofs = [100, 500, 1000, 2000] << u.s
    # Propagate the secondary body to generate its position coordinates.
    rr, vv = cowell(
        Earth.k,
        orbit.r,
        orbit.v,
        tofs,
    )
    pos_coords = rr  # Trajectory of the secondary body.

    orb = Orbit.from_classical(
        attractor=Earth,
        a=16000 * u.km,
        ecc=0.53 * u.one,
        inc=5 * u.deg,
        raan=5 * u.deg,
        argp=10 * u.deg,
        nu=30 * u.deg
    )

    los_event = LosEvent(Earth, pos_coords, terminal=True)
    events = [los_event]
    tofs = [1, 5, 10, 100, 1000, 2000, 3000, 5000] << u.s

    r, v = cowell(
        Earth.k,
        orb.r,
        orb.v,
        tofs,
        events=events,
    )

    assert_quantity_allclose(los_event.last_t, t_los)
