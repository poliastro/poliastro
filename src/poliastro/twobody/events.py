from astropy import units as u
from numpy.linalg import norm


class LithobrakeEvent:
    def __init__(self, R):
        self._R = R
        self._last_t = None

    @property
    def terminal(self):
        # tell scipy to stop the integration at H = R (impact)
        return True

    @property
    def last_t(self):
        return self._last_t * u.s

    def __call__(self, t, u):
        self._last_t = t
        H = norm(u[:3])
        # scipy will search for H - R = 0
        return H - self._R
