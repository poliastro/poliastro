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
    """Converts from classical orbital elements to vectors.

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
    # FIXME: rv2coe and coe2rv are not transitive
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


def coe2mee(p, ecc, inc, raan, argp, nu):
    """Converts from classical orbital elements to modified equinoctial
    orbital elements.

    The definition of the modified equinoctial orbital elements is taken from
    [Walker, 1985].

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

    Notes
    -----
    The conversion equations are taken directly from the original paper.

    """
    lonper = raan + argp
    f = ecc * np.cos(lonper)
    g = ecc * np.sin(lonper)
    # TODO: Check polar case (see [Walker, 1985])
    h = np.tan(inc / 2) * np.cos(raan)
    k = np.tan(inc / 2) * np.sin(raan)
    L = lonper + nu
    return p, f, g, h, k, L


def mee2coe(p, f, g, h, k, L):
    """Converts from modified equinoctial orbital elements to classical
    orbital elements.

    The definition of the modified equinoctial orbital elements is taken from
    [Walker, 1985].

    Notes
    -----
    The conversion is always safe because arctan2 works also for 0, 0
    arguments.

    """
    ecc = np.sqrt(f**2 + g**2)
    inc = 2 * np.arctan(np.sqrt(h**2 + k**2))
    lonper = np.arctan2(g, f)
    raan = np.arctan2(k, h) % (2 * np.pi)
    argp = (lonper - raan) % (2 * np.pi)
    nu = (L - lonper) % (2 * np.pi)
    return p, ecc, inc, raan, argp, nu
