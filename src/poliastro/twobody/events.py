import numpy as np
from astropy import units as u
from astropy.coordinates import get_body_barycentric_posvel
from numpy.linalg import norm

from poliastro.core.events import (
    eclipse_function as eclipse_function_fast,
    visibility_function as visibility_function_fast,
)
from poliastro.core.spheroid_location import (
    cartesian_to_ellipsoidal as cartesian_to_ellipsoidal_fast,
)
from poliastro.spheroid_location import SpheroidLocation


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
    terminal: bool, optional
        Whether to terminate integration if this event occurs, defaults to True.
    direction: float, optional
        Handle triggering of event based on whether latitude is crossed from above
        or below, defaults to 0, i.e., event is triggered while traversing from both directions.

    """

    def __init__(self, orbit, lat, terminal=False, direction=0):
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
        self.k = self._primary_body.k.to_value(u.km ** 3 / u.s ** 2)
        self.R_sec = self._secondary_body.R.to_value(u.km)
        self.R_primary = self._primary_body.R.to_value(u.km)

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

        r_sec = super().__call__(t, u_, k)
        shadow_function = eclipse_function_fast(
            self.k,
            u_,
            r_sec,
            self.R_sec,
            self.R_primary,
            umbra=False,
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

        r_sec = super().__call__(t, u_, k)
        shadow_function = eclipse_function_fast(
            self.k, u_, r_sec, self.R_sec, self.R_primary
        )

        return shadow_function


class NodeCrossEvent(Event):
    """Detect equatorial node (ascending or descending) crossings.

    Parameters
    ----------
    terminal: bool, optional
        Whether to terminate integration when the event occurs, defaults to False.
    direction: float, optional
        Handle triggering of event based on whether the node is crossed from above
        i.e. descending node, or is crossed from below i.e. ascending node, defaults to 0,
        i.e. event is triggered during both crossings.

    """

    def __init__(self, terminal=False, direction=0):
        super().__init__(terminal, direction)

    def __call__(self, t, u_, k):
        self._last_t = t
        # Check if the z coordinate of the satellite is zero.
        return u_[2]


class SatelliteVisibility(Event):
    """Detects whether a satellite is visible from a ground station.

    Parameters
    ----------
    el: float
        Minimum elevation angle of the station.
    lon: float
        Longitude of the station.
    lat: float
        Latitude of the station.
    h: float
        Geodetic height.
    body: poliastro.bodies.Body
        The body on which the station is defined.

    """

    def __init__(self, el, lon, lat, h, body, terminal=False, direction=0):
        super().__init__(terminal, direction)
        self._el = el.to(u.rad).value
        sph_loc = SpheroidLocation(lon, lat, h, body)
        self._lon = lon.to(u.rad).value  # Unused
        self._lat = lat.to(u.rad).value
        self._H = h.to(u.km).value
        self._R = body.R.to(u.km).value
        self._N = (
            sph_loc.N << u.km
        ).value  # Instead compute in core to prevent creating SpheroidLocation object?
        self._k = body.k.to_value(u.km ** 3 / u.s ** 2)

    def __call__(self, t, u_, k):
        self._last_t = t

        visibility_function = visibility_function_fast(
            self._k,
            u_,
            self._N,
            self._lat,
            self._H,
            self._R,
            self._el,
        )
        return visibility_function
