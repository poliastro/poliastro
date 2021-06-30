import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import (
    GCRS,
    ITRS,
    CartesianRepresentation,
    SphericalRepresentation,
)
from astropy.tests.helper import assert_quantity_allclose
from numpy.linalg import norm

from poliastro.bodies import Earth
from poliastro.constants import H0_earth, rho0_earth
from poliastro.core.perturbations import atmospheric_drag_exponential
from poliastro.core.propagation import func_twobody
from poliastro.twobody import Orbit
from poliastro.twobody.events import (
    AltitudeCrossEvent,
    LatitudeCrossEvent,
    LongitudeCrossEvent,
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


def test_latitude_longitude_cross_event():
    R = Earth.R.to(u.km).value
    orbit = Orbit.circular(Earth, 230 * u.km)

    t_lat = 1989.9932 * u.s

    thresh_lat = 0 * u.deg
    latitude_cross_event = LatitudeCrossEvent(thresh_lat.value, R)

    tofs = [5] * u.d

    events = [latitude_cross_event]
    rr, _ = cowell(
        Earth.k,
        orbit.r,
        orbit.v,
        tofs,
        events=events,
    )

    obstime = orbit.epoch + t_lat
    gcrs_xyz = GCRS(
        rr[-1],
        obstime=obstime,
        representation_type=CartesianRepresentation,
    )
    itrs_xyz = gcrs_xyz.transform_to(ITRS(obstime=obstime))
    itrs_latlon_pos = itrs_xyz.represent_as(SphericalRepresentation)
    orbit_lat = itrs_latlon_pos.lat.to(u.deg)

    assert_quantity_allclose(latitude_cross_event.last_t, t_lat)
    assert_quantity_allclose(orbit_lat, thresh_lat, atol=1.2e-4 * u.deg)


def test_longitude_cross():
    R = Earth.R.to(u.km).value
    orbit = Orbit.circular(Earth, 230 * u.km)

    t_lon = 0.06289492 * u.s

    thresh_lon = 79.810036 * u.deg
    longitude_cross_event = LongitudeCrossEvent(thresh_lon.value, R)

    tofs = [5] * u.d

    events = [longitude_cross_event]
    rr, _ = cowell(
        Earth.k,
        orbit.r,
        orbit.v,
        tofs,
        events=events,
    )

    obstime = orbit.epoch + t_lon
    gcrs_xyz = GCRS(
        rr[-1],
        obstime=obstime,
        representation_type=CartesianRepresentation,
    )
    itrs_xyz = gcrs_xyz.transform_to(ITRS(obstime=obstime))
    itrs_latlon_pos = itrs_xyz.represent_as(SphericalRepresentation)
    orbit_lon = itrs_latlon_pos.lon.to(u.deg)

    assert_quantity_allclose(longitude_cross_event.last_t, t_lon)
    assert_quantity_allclose(orbit_lon, thresh_lon)
