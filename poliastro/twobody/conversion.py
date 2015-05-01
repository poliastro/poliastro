# coding: utf-8
"""Conversion functions.

"""

import numpy as np
from numpy import cos, sin, sqrt

from astropy import units as u

from poliastro.util import transform, norm


def rv_pqw(k, p, ecc, nu):
    """Returns r and v vectors in perifocal frame.

    """
    r_pqw = (np.array([cos(nu), sin(nu), 0 * nu]) * p / (1 + ecc * cos(nu))).T
    v_pqw = (np.array([-sin(nu), (ecc + cos(nu)), 0]) * sqrt(k / p)).T
    return r_pqw, v_pqw


def coe2rv(k, p, ecc, inc, raan, argp, nu):
    """Converts from orbital elements to vectors.

    Parameters
    ----------
    k : float
        Standard gravitational parameter (km^3 / s^2).
    p : float
        Semi-latus rectum or parameter (km).
    ecc : float
        Eccentricity.
    inc : float
        Inclination (rad).
    omega : float
        Longitude of ascending node (rad).
    argp : float
        Argument of perigee (rad).
    nu : float
        True anomaly (rad).

    """
    r_pqw, v_pqw = rv_pqw(k, p, ecc, nu)

    r_ijk = transform(r_pqw, -argp, 'z', u.rad)
    r_ijk = transform(r_ijk, -inc, 'x', u.rad)
    r_ijk = transform(r_ijk, -raan, 'z', u.rad)
    v_ijk = transform(v_pqw, -argp, 'z', u.rad)
    v_ijk = transform(v_ijk, -inc, 'x', u.rad)
    v_ijk = transform(v_ijk, -raan, 'z', u.rad)

    return r_ijk, v_ijk


def rv2coe(k, r, v, tol=1e-8):
    """Converts from vectors to orbital elements.

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
    a = p / (1 - ecc ** 2)
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

    return a, ecc, inc, raan, argp, nu
