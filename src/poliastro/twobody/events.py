from astropy import units as u
from numpy.linalg import norm


class LithobrakeEvent:
    """Terminal event that detects impact with the attractor surface.

    Parameters
    ----------
    R : float
        Radius of the attractor.

    """

    def __init__(self, R):
        self._R = R
        self._last_t = None

    @property
    def terminal(self):
        # Tell SciPy to stop the integration at H = R (impact)
        return True

    @property
    def last_t(self):
        return self._last_t * u.s

    def __call__(self, t, u, k):
        self._last_t = t
        H = norm(u[:3])
        # SciPy will search for H - R = 0
        print(H - self._R)
        return H - self._R


class AltitudeCrossEvent:
    """Detect if a satellite crosses a specific threshold altitude.

    Parameters
    ----------
    R: ~astropy.units.Quantity
        Radius of the attractor (km).
    thresh_H: ~astropy.units.Quantity
        Threshold altitude (in km), defaults to 100 km.
    terminal: bool
        Whether to terminate integration if this event occurs, defaults to True.

    """
    def __init__(self, R, thresh_H=100*u.km, terminal=True):
        self._R = R.to(u.km).value
        self._thresh_H = thresh_H.to(u.km).value  # Threshold height from the ground.
        self._terminal = terminal
        self._last_t = None

    @property
    def terminal(self):
        # Orekit's API stops propagation when descending, but not when ascending.
        return self._terminal

    @property
    def last_t(self):
        return self._last_t * u.s

    def __call__(self, t, u, k):
        self._last_t = t
        H = norm(u[:3])
        # H is from the center of the attractor.
        return H - self._R - self._thresh_H  # If this goes from +ve to -ve, altitude is decreasing.
