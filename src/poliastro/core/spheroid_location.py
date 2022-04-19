""" Low level calculations for oblate spheroid locations """

import numpy as np
from numba import njit as jit

from poliastro._math.linalg import norm


@jit
def cartesian_cords(_a, _c, _lon, _lat, _h):
    """Calculates cartesian coordinates.

    Parameters
    ----------
    _a : float
        Semi-major axis
    _c : float
        Semi-minor axis
    _lon : float
        Geodetic longitude
    _lat : float
        Geodetic latitude
    _h : float
        Geodetic height

    """
    e2 = 1 - (_c / _a) ** 2
    N = _a / np.sqrt(1 - e2 * np.sin(_lon) ** 2)

    x = (N + _h) * np.cos(_lon) * np.cos(_lat)
    y = (N + _h) * np.cos(_lon) * np.sin(_lat)
    z = ((1 - e2) * N + _h) * np.sin(_lon)
    return x, y, z


@jit
def f(_a, _c):
    """Get first flattening.

    Parameters
    ----------
    _a : float
        Semi-major axis
    _c : float
        Semi-minor axis

    """
    return 1 - _c / _a


@jit
def N(a, b, c, cartesian_cords):
    """Normal vector of the ellipsoid at the given location.

    Parameters
    ----------
    a : float
        Semi-major axis
    b : float
        Equatorial radius
    c : float
        Semi-minor axis
    cartesian_cords : numpy.ndarray
        Cartesian coordinates

    """
    x, y, z = cartesian_cords
    N = np.array([2 * x / a**2, 2 * y / b**2, 2 * z / c**2])
    N /= norm(N)
    return N


@jit
def tangential_vecs(N):
    """Returns orthonormal vectors tangential to the ellipsoid at the given location.

    Parameters
    ----------
    N : numpy.ndarray
        Normal vector of the ellipsoid

    """
    u = np.array([1.0, 0, 0])
    u -= (u @ N) * N
    u /= norm(u)
    v = np.cross(N, u)

    return u, v


@jit
def radius_of_curvature(_a, _c, _lat):
    """Radius of curvature of the meridian at the latitude of the given location.

    Parameters
    ----------
    _a : float
        Semi-major axis
    _c : float
        Semi-minor axis
    _lat : float
        Geodetic latitude

    """
    e2 = 1 - (_c / _a) ** 2
    rc = _a * (1 - e2) / (1 - e2 * np.sin(_lat) ** 2) ** 1.5
    return rc


@jit
def distance(cartesian_cords, px, py, pz):
    """Calculates the distance from an arbitrary point to the given location (Cartesian coordinates).

    Parameters
    ----------
    cartesian_cords : numpy.ndarray
        Cartesian coordinates
    px : float
        x-coordinate of the point
    py : float
        y-coordinate of the point
    pz : float
        z-coordinate of the point

    """
    c = cartesian_cords
    u = np.array([px, py, pz])
    d = norm(c - u)
    return d


@jit
def is_visible(cartesian_cords, px, py, pz, N):
    """Determine whether an object located at a given point is visible from the given location.

    Parameters
    ----------
    cartesian_cords : numpy.ndarray
        Cartesian coordinates
    px : float
        x-coordinate of the point
    py : float
        y-coordinate of the point
    pz : float
        z-coordinate of the point
    N : numpy.ndarray
        Normal vector of the ellipsoid at the given location.

    """
    c = cartesian_cords
    u = np.array([px, py, pz])

    d = -(N @ c)
    p = (N @ u) + d
    return p >= 0


@jit
def cartesian_to_ellipsoidal(_a, _c, x, y, z):
    """
    Converts cartesian coordinates to ellipsoidal coordinates for the given ellipsoid.
    Instead of the iterative formula, the function uses the approximation introduced in
    Bowring, B. R. (1976). TRANSFORMATION FROM SPATIAL TO GEOGRAPHICAL COORDINATES

    Parameters
    ----------
    _a : float
        Semi-major axis
    _c : float
        Semi-minor axis
    x : float
        x coordinate
    y : float
        y coordinate
    z : float
        z coordinate

    """
    e2 = 1 - (_c / _a) ** 2
    e2_ = e2 / (1 - e2)
    p = np.sqrt(x**2 + y**2)
    th = np.arctan(z * _a / (p * _c))
    lon = np.arctan2(
        y, x
    )  # Use `arctan2` so that lon lies in the range: [-pi, +pi]
    lat = np.arctan(
        (z + e2_ * _c * np.sin(th) ** 3) / (p - e2 * _a * np.cos(th) ** 3)
    )

    v = _a / np.sqrt(1 - e2 * np.sin(lat) ** 2)
    h = x / np.cos(lat) - v if lat == 0.0 else z / np.sin(lat) - (1 - e2) * v

    return lat, lon, h
