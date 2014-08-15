# coding: utf-8
"""Conversion functions.

"""

import numpy as np
from numpy import cos, sin, sqrt

from astropy import units as u
u.one = u.dimensionless_unscaled  # astropy #1980

from poliastro.util import transform
from poliastro.util import norm


def rv_pqw(k, p, ecc, nu):
    """Returns r and v vectors in perifocal frame.

    """
    r_pqw = (np.array([cos(nu), sin(nu), 0 * nu]) * p / (1 + ecc * cos(nu))).T
    v_pqw = (np.array([-sin(nu), (ecc + cos(nu)), 0]) * sqrt(k / p)).T
    return r_pqw, v_pqw


def coe2rv(k, a, ecc, inc, raan, argp, nu):
    """Converts from orbital elements to vectors.

    Parameters
    ----------
    k : float
        Standard gravitational parameter (km^3 / s^2).
    a : float
        Semi-major axis (km).
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
    p = a * (1 - ecc ** 2)
    r_pqw, v_pqw = rv_pqw(k, p, ecc, nu)

    r_ijk = transform(r_pqw, -argp, 'z', u.rad)
    r_ijk = transform(r_ijk, -inc, 'x', u.rad)
    r_ijk = transform(r_ijk, -raan, 'z', u.rad)
    v_ijk = transform(v_pqw, -argp, 'z', u.rad)
    v_ijk = transform(v_ijk, -inc, 'x', u.rad)
    v_ijk = transform(v_ijk, -raan, 'z', u.rad)

    return r_ijk, v_ijk


def rv2coe(k, r, v):
    """Converts from vectors to orbital elements.

    Parameters
    ----------
    k : float
        Standard gravitational parameter (km^3 / s^2).
    r : array
        Position vector (km).
    v : array
        Velocity vector (km / s).

    """
    h = np.cross(r, v)
    n = np.cross([0, 0, 1], h) / norm(h)
    e = ((v.dot(v) - k / (norm(r))) * r - r.dot(v) * v) / k
    ecc = norm(e)
    p = h.dot(h) / k
    a = p / (1 - ecc ** 2)

    inc = np.arccos(h[2] / norm(h))
    raan = np.arctan2(n[1], n[0]) % (2 * np.pi)
    argp = np.arctan2(h.dot(np.cross(n, e)) / norm(h), e.dot(n)) % (2 * np.pi)
    nu = np.arctan2(h.dot(np.cross(e, r)) / norm(h), r.dot(e)) % (2 * np.pi)

    return a, ecc, inc, raan, argp, nu
