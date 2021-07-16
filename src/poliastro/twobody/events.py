import numpy as np
from astropy import units as u
from numpy.linalg import norm

from poliastro.core.spheroid_location import (
    cartesian_to_ellipsoidal as cartesian_to_ellipsoidal_fast,
)


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
        Threshold latitude (in degrees).
    terminal: bool
        Whether to terminate integration if this event occurs, defaults to True.
    direction: float
        Handle triggering of event based on whether latitude is crossed from above
        or below, defaults to 0, i.e., event is triggered while traversing from both directions.

    """

    def __init__(self, orbit, lat, terminal=True, direction=0):
        super().__init__(terminal, direction)

        self._R = orbit.attractor.R
        self._R_polar = orbit.attractor.R_polar
        self._epoch = orbit.epoch
        self._lat = lat.to(u.deg).value  # Threshold latitude (in degrees).

    def __call__(self, t, u_, k):
        self._last_t = t
        pos_on_body = (u_[:3] / norm(u_[:3])) * self._R
        lat_, _, _ = cartesian_to_ellipsoidal_fast(self._R, self._R_polar, *pos_on_body)

        return np.rad2deg(lat_) - self._lat
