# coding: utf-8
"""Functions to define orbits from position and velocity vectors.

"""
import numpy as np
from astropy import units as u

from poliastro.util import norm

from poliastro.twobody.core import State


def rv2coe(k, r, v, tol=1e-8):
    """Converts from vectors to classical orbital elements.

    Parameters
    ----------
    k : float
        Standard gravitational parameter (km^3 / s^2).
    r : array
        Position vector (km).
    v : array
        Velocity vector (km / s).
    tol : float, optional
        Tolerance for eccentricity and inclination checks, default to 1e-8.

    """
    h = np.cross(r, v)
    n = np.cross([0, 0, 1], h) / norm(h)
    e = ((v.dot(v) - k / (norm(r))) * r - r.dot(v) * v) / k
    ecc = norm(e)
    p = h.dot(h) / k
    inc = np.arccos(h[2] / norm(h))

    circular = ecc < tol
    equatorial = abs(inc) < tol

    if equatorial and not circular:
        raan = 0
        argp = np.arctan2(e[1], e[0]) % (2 * np.pi)  # Longitude of periapsis
        nu = (np.arctan2(h.dot(np.cross(e, r)) / norm(h), r.dot(e)) %
              (2 * np.pi))
    elif not equatorial and circular:
        raan = np.arctan2(n[1], n[0]) % (2 * np.pi)
        argp = 0
        # Argument of latitude
        nu = (np.arctan2(r.dot(np.cross(h, n)) / norm(h), r.dot(n)) %
              (2 * np.pi))
    elif equatorial and circular:
        raan = 0
        argp = 0
        nu = np.arctan2(r[1], r[0]) % (2 * np.pi)  # True longitude
    else:
        raan = np.arctan2(n[1], n[0]) % (2 * np.pi)
        argp = (np.arctan2(e.dot(np.cross(h, n)) / norm(h), e.dot(n)) %
                (2 * np.pi))
        nu = (np.arctan2(r.dot(np.cross(h, e)) / norm(h), r.dot(e))
              % (2 * np.pi))

    return p, ecc, inc, raan, argp, nu


class RVState(State):
    def __init__(self, attractor, r, v, epoch):
        super(RVState, self).__init__(attractor, epoch)
        self._r = r
        self._v = v

    @property
    def r(self):
        return self._r

    @property
    def v(self):
        return self._v

    def to_vectors(self):
        return self

    def to_classical(self):
        (p, ecc, inc, raan, argp, nu
         ) = rv2coe(self.attractor.k.to(u.km ** 3 / u.s ** 2).value,
                    self.r.to(u.km).value,
                    self.v.to(u.km / u.s).value)

        return super(RVState, self).from_classical(self.attractor,
                                                   p * u.km,
                                                   ecc * u.one,
                                                   inc * u.rad,
                                                   raan * u.rad,
                                                   argp * u.rad,
                                                   nu * u.rad,
                                                   self.epoch)
