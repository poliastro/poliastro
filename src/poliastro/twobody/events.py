import numpy as np
from astropy import units as u
from astropy.coordinates import get_body_barycentric_posvel
from numpy.linalg import norm

from poliastro.core.elements import rv2coe
from poliastro.core.spheroid_location import (
    cartesian_to_ellipsoidal as cartesian_to_ellipsoidal_fast,
)
from poliastro.twobody.orbit import Orbit


class Event:
    """Base class for event functionalities.

    Parameters
    ----------
    terminal: bool
        Whether to terminate integration if this event occurs.
    direction: float
        Handle triggering of event.

    """

    def __init__(self, terminal, direction):
        self._terminal, self._direction = terminal, direction
        self._last_t = None

    @property
    def terminal(self):
        return self._terminal

    @property
    def direction(self):
        return self._direction

    @property
    def last_t(self):
        return self._last_t * u.s

    def __call__(self, t, u, k):
        raise NotImplementedError


class AltitudeCrossEvent(Event):
    """Detect if a satellite crosses a specific threshold altitude.

    Parameters
    ----------
    alt: float
        Threshold altitude (km).
    R: float
        Radius of the attractor (km).
    terminal: bool
        Whether to terminate integration if this event occurs.
    direction: float
        Handle triggering of event based on whether altitude is crossed from above
        or below, defaults to -1, i.e., event is triggered only if altitude is
        crossed from above (decreasing altitude).

    """

    def __init__(self, alt, R, terminal=True, direction=-1):
        super().__init__(terminal, direction)
        self._R = R
        self._alt = alt  # Threshold altitude from the ground.

    def __call__(self, t, u, k):
        self._last_t = t
        r_norm = norm(u[:3])

        return (
            r_norm - self._R - self._alt
        )  # If this goes from +ve to -ve, altitude is decreasing.


class LithobrakeEvent(AltitudeCrossEvent):
    """Terminal event that detects impact with the attractor surface.

    Parameters
    ----------
    R : float
        Radius of the attractor (km).
    terminal: bool
        Whether to terminate integration if this event occurs.

    """

    def __init__(self, R, terminal=True):
        super().__init__(0, R, terminal, direction=-1)


class LatitudeCrossEvent(Event):
    """Detect if a satellite crosses a specific threshold latitude.

    Parameters
    ----------
    orbit: ~poliastro.twobody.orbit.Orbit
        Orbit.
    lat: astropy.quantity.Quantity
        Threshold latitude.
    terminal: bool
        Whether to terminate integration if this event occurs, defaults to True.
    direction: float
        Handle triggering of event based on whether latitude is crossed from above
        or below, defaults to 0, i.e., event is triggered while traversing from both directions.

    """

    def __init__(self, orbit, lat, terminal=True, direction=0):
        super().__init__(terminal, direction)

        self._R = orbit.attractor.R.to(u.m).value
        self._R_polar = orbit.attractor.R_polar.to(u.m).value
        self._epoch = orbit.epoch
        self._lat = lat.to(u.deg).value  # Threshold latitude (in degrees).

    def __call__(self, t, u_, k):
        self._last_t = t
        pos_on_body = (u_[:3] / norm(u_[:3])) * self._R
        lat_, _, _ = cartesian_to_ellipsoidal_fast(self._R, self._R_polar, *pos_on_body)

        return np.rad2deg(lat_) - self._lat


class EclipseEvent(Event):
    """Base class for the eclipse event.

    Parameters
    ----------
    orbit: poliastro.twobody.orbit.Orbit
        Orbit of the satellite.
    terminal: bool, optional
        Whether to terminate integration when the event occurs, defaults to False.
    direction: float, optional
        Specify which direction must the event trigger, defaults to 0.

    """

    def __init__(self, orbit, terminal=False, direction=0):
        super().__init__(terminal, direction)
        self._primary_body = orbit.attractor
        self._secondary_body = orbit.attractor.parent
        self._epoch = orbit.epoch

    def __call__(self, t, u_, k):
        # Solve for primary and secondary bodies position w.r.t. solar system
        # barycenter at a particular epoch.
        (r_primary_wrt_ssb, _), (r_secondary_wrt_ssb, _) = [
            get_body_barycentric_posvel(body.name, self._epoch + t * u.s)
            for body in (self._primary_body, self._secondary_body)
        ]
        r_sec = ((r_secondary_wrt_ssb - r_primary_wrt_ssb).xyz << u.km).value

        return r_sec


class PenumbraEvent(EclipseEvent):
    """Detect whether a satellite is in penumbra or not.

    Parameters
    ----------
    orbit: poliastro.twobody.orbit.Orbit
        Orbit of the satellite.
    terminal: bool, optional
        Whether to terminate integration when the event occurs, defaults to False.
    direction: float, optional
        Handle triggering of event based on whether entry is into or out of
        penumbra, defaults to 0, i.e., event is triggered at both, entry and exit points.

    """

    def __init__(self, orbit, terminal=False, direction=0):
        super().__init__(orbit, terminal, direction)

    def __call__(self, t, u_, k):
        self._last_t = t

        k = self._primary_body.k.to_value(u.km ** 3 / u.s ** 2)
        R_sec = self._secondary_body.R.to_value(u.km)
        R_primary = self._primary_body.R.to_value(u.km)

        r_sec = super().__call__(t, u_, k)

        p, ecc, inc, raan, argp, nu = rv2coe(k, u_[:3], u_[3:])

        # A dummy orbit is needed to access the `pqw` method.
        dummy_orbit = Orbit.from_vectors(
            self._primary_body, u_[:3] << u.km, u_[3:] << u.km / u.s
        )
        P_, Q_, _ = dummy_orbit.pqw()

        # Px = np.cos(raan) * np.cos(argp) - np.sin(raan) * np.sin(argp) * np.cos(inc)
        # Py = np.cos(raan) * np.sin(argp) + np.sin(raan) * np.cos(argp) * np.cos(inc)
        # Pz = np.sin(raan) * np.sin(inc)
        # P_ = np.array([Px, Py, Pz])
        # Qx = -np.sin(raan) * np.cos(argp) - np.cos(raan) * np.sin(argp) * np.cos(inc)
        # Qy = -np.sin(raan) * np.sin(argp) + np.cos(raan) * np.cos(argp) * np.cos(inc)
        # Qz = np.cos(raan) * np.sin(inc)
        # Q_ = np.array([Qx, Qy, Qz])

        r_sec_norm = norm(r_sec)
        beta = np.dot(P_, r_sec) / r_sec_norm
        zeta = np.dot(Q_, r_sec) / r_sec_norm

        delta_p = np.arcsin((R_sec + R_primary) / r_sec_norm)

        cos_psi = beta * np.cos(nu) + zeta * np.sin(nu)
        shadow_function = (
            ((R_primary ** 2) * (1 + ecc * np.cos(nu)) ** 2)
            + (p ** 2) * (cos_psi ** 2)
            - p ** 2
            - (2 * p * R_primary * cos_psi) * (1 + ecc * np.cos(nu)) * np.sin(delta_p)
        )

        return shadow_function


class UmbraEvent(EclipseEvent):
    """Detect whether a satellite is in umbra or not.

    Parameters
    ----------
    orbit: poliastro.twobody.orbit.Orbit
        Orbit of the satellite.
    terminal: bool, optional
        Whether to terminate integration when the event occurs, defaults to False.
    direction: float, optional
        Handle triggering of event based on whether entry is into or out of
        umbra, defaults to 0, i.e., event is triggered at both, entry and exit points.

    """

    def __init__(self, orbit, terminal=False, direction=0):
        super().__init__(orbit, terminal, direction)

    def __call__(self, t, u_, k):
        self._last_t = t

        k = self._primary_body.k.to_value(u.km ** 3 / u.s ** 2)
        R_sec = self._secondary_body.R.to_value(u.km)
        R_primary = self._primary_body.R.to_value(u.km)

        r_sec = super().__call__(t, u_, k)

        p, ecc, inc, raan, argp, nu = rv2coe(k, u_[:3], u_[3:])

        dummy_orbit = Orbit.from_vectors(
            self._primary_body, u_[:3] << u.km, u_[3:] << u.km / u.s
        )
        P_, Q_, _ = dummy_orbit.pqw()

        # Px = np.cos(raan) * np.cos(argp) - np.sin(raan) * np.sin(argp) * np.cos(inc)
        # Py = np.cos(raan) * np.sin(argp) + np.sin(raan) * np.cos(argp) * np.cos(inc)
        # Pz = np.sin(raan) * np.sin(inc)
        # P_ = np.array([Px, Py, Pz])
        # Qx = -np.sin(raan) * np.cos(argp) - np.cos(raan) * np.sin(argp) * np.cos(inc)
        # Qy = -np.sin(raan) * np.sin(argp) + np.cos(raan) * np.cos(argp) * np.cos(inc)
        # Qz = np.cos(raan) * np.sin(inc)
        # Q_ = np.array([Qx, Qy, Qz])

        r_sec_norm = norm(r_sec)
        beta = np.dot(P_, r_sec) / r_sec_norm
        zeta = np.dot(Q_, r_sec) / r_sec_norm

        delta_u = np.arcsin((R_sec - R_primary) / r_sec_norm)

        cos_psi = beta * np.cos(nu) + zeta * np.sin(nu)
        shadow_function = (
            ((R_primary ** 2) * (1 + ecc * np.cos(nu)) ** 2)
            + (p ** 2) * (cos_psi ** 2)
            - p ** 2
            + (2 * p * R_primary * cos_psi) * (1 + ecc * np.cos(nu)) * np.sin(delta_u)
        )

        return shadow_function
