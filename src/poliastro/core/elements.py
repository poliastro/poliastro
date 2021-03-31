"""This module contains a set of functions that can be used to
convert between different elements that define the orbit of a body.
"""

import numpy as np
from numba import njit as jit, prange
from numpy import cos, cross, sin, sqrt
from numpy.linalg import norm

from .angles import E_to_nu, F_to_nu
from .util import rotation_matrix


@jit
def rv_pqw(k, p, ecc, nu):
    r"""Returns r and v vectors in perifocal frame.

    Parameters
    ----------
    k : float
        Standard gravitational parameter (km^3 / s^2).
    p : float
        Semi-latus rectum or parameter (km).
    ecc : float
        Eccentricity.
    nu: float
        True anomaly (rad).

    Returns
    -------

    r: ndarray
        Position. Dimension 3 vector
    v: ndarray
        Velocity. Dimension 3 vector

    Notes
    -----
    These formulas can be checked at Curtis 3rd. Edition, page 110. Also the
    example proposed is 2.11 of Curtis 3rd Edition book.

    .. math::

        \vec{r} = \frac{h^2}{\mu}\frac{1}{1 + e\cos(\theta)}\begin{bmatrix}
        \cos(\theta)\\
        \sin(\theta)\\
        0
        \end{bmatrix} \\\\\\

        \vec{v} = \frac{h^2}{\mu}\begin{bmatrix}
        -\sin(\theta)\\
        e+\cos(\theta)\\
        0
        \end{bmatrix}

    Examples
    --------
    >>> from poliastro.constants import GM_earth
    >>> k = GM_earth.value  # Earth gravitational parameter
    >>> ecc = 0.3  # Eccentricity
    >>> h = 60000e6  # Angular momentum of the orbit (m**2 / s)
    >>> nu = np.deg2rad(120)  # True Anomaly (rad)
    >>> p = h**2 / k  # Parameter of the orbit
    >>> r, v = rv_pqw(k, p, ecc, nu)
    >>> # Printing the results
    r = [-5312706.25105345  9201877.15251336    0] [m]
    v = [-5753.30180931 -1328.66813933  0] [m]/[s]

    """
    pqw = np.array([[cos(nu), sin(nu), 0], [-sin(nu), ecc + cos(nu), 0]]) * np.array(
        [[p / (1 + ecc * cos(nu))], [sqrt(k / p)]]
    )
    return pqw


@jit
def coe_rotation_matrix(inc, raan, argp):
    """Create a rotation matrix for coe transformation"""
    r = rotation_matrix(raan, 2)
    r = r @ rotation_matrix(inc, 0)
    r = r @ rotation_matrix(argp, 2)
    return r


@jit
def coe2rv(k, p, ecc, inc, raan, argp, nu):
    r"""Converts from classical orbital to state vectors.

    Classical orbital elements are converted into position and velocity
    vectors by `rv_pqw` algorithm. A rotation matrix is applied to position
    and velocity vectors to get them expressed in terms of an IJK basis.

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

    Returns
    -------
    r_ijk: np.array
        Position vector in basis ijk.
    v_ijk: np.array
        Velocity vector in basis ijk.

    Notes
    -----

    .. math::
        \begin{align}
            \vec{r}_{IJK} &= [ROT3(-\Omega)][ROT1(-i)][ROT3(-\omega)]\vec{r}_{PQW}
                               = \left [ \frac{IJK}{PQW} \right ]\vec{r}_{PQW}\\
            \vec{v}_{IJK} &= [ROT3(-\Omega)][ROT1(-i)][ROT3(-\omega)]\vec{v}_{PQW}
                               = \left [ \frac{IJK}{PQW} \right ]\vec{v}_{PQW}\\
        \end{align}

    Previous rotations (3-1-3) can be expressed in terms of a single rotation matrix:

    .. math::
        \left [ \frac{IJK}{PQW} \right ]

    .. math::
        \begin{bmatrix}
        \cos(\Omega)\cos(\omega) - \sin(\Omega)\sin(\omega)\cos(i) & -\cos(\Omega)\sin(\omega) - \sin(\Omega)\cos(\omega)\cos(i) & \sin(\Omega)\sin(i)\\
        \sin(\Omega)\cos(\omega) + \cos(\Omega)\sin(\omega)\cos(i) & -\sin(\Omega)\sin(\omega) + \cos(\Omega)\cos(\omega)\cos(i) & -\cos(\Omega)\sin(i)\\
        \sin(\omega)\sin(i) & \cos(\omega)\sin(i) & \cos(i)
        \end{bmatrix}

    """
    pqw = rv_pqw(k, p, ecc, nu)
    rm = coe_rotation_matrix(inc, raan, argp)

    ijk = pqw @ rm.T

    return ijk


@jit(parallel=True)
def coe2rv_many(k, p, ecc, inc, raan, argp, nu):

    n = nu.shape[0]
    rr = np.zeros((n, 3))
    vv = np.zeros((n, 3))

    for i in prange(n):
        rr[i, :], vv[i, :] = coe2rv(k[i], p[i], ecc[i], inc[i], raan[i], argp[i], nu[i])

    return rr, vv


@jit
def coe2mee(p, ecc, inc, raan, argp, nu):
    r"""Converts from classical orbital elements to modified equinoctial
    orbital elements.

    The definition of the modified equinoctial orbital elements is taken from
    [Walker, 1985].

    The modified equinoctial orbital elements are a set of orbital elements that are useful for
    trajectory analysis and optimization. They are valid for circular, elliptic, and hyperbolic
    orbits. These direct modified equinoctial equations exhibit no singularity for zero
    eccentricity and orbital inclinations equal to 0 and 90 degrees. However, two of the
    components are singular for an orbital inclination of 180 degrees.

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

    Returns
    -------
    p: float
        Semi-latus rectum or parameter
    f: float
        Equinoctial parameter f
    g: float
        Equinoctial parameter g
    h: float
        Equinoctial parameter h
    k: float
        Equinoctial parameter k
    L: float
        Longitude

    Notes
    -----
    The conversion equations are taken directly from the original paper:

    .. math::
        \begin{align}
        p &= a(1-e^2) \\
        f &= e\cos(\omega + \Omega) \\
        g &= e\sin(\omega + \Omega) \\
        h &= \tan(\frac{i}{2})\cos(\Omega) \\
        k &= \tan(\frac{i}{2})\sin(\Omega) \\
        L &= \Omega + \omega + \theta \\
        \end{align}

    """
    lonper = raan + argp
    f = ecc * np.cos(lonper)
    g = ecc * np.sin(lonper)
    # TODO: Check polar case (see [Walker, 1985])
    h = np.tan(inc / 2) * np.cos(raan)
    k = np.tan(inc / 2) * np.sin(raan)
    L = lonper + nu
    return p, f, g, h, k, L


@jit
def rv2coe(k, r, v, tol=1e-8):
    r"""Converts from vectors to classical orbital elements.

    Parameters
    ----------
    k : float
        Standard gravitational parameter (km^3 / s^2)
    r : array
        Position vector (km)
    v : array
        Velocity vector (km / s)
    tol : float, optional
        Tolerance for eccentricity and inclination checks, default to 1e-8

    Returns
    -------
    p : float
        Semi-latus rectum of parameter (km)
    ecc: float
        Eccentricity
    inc: float
        Inclination (rad)
    raan: float
        Right ascension of the ascending nod (rad)
    argp: float
        Argument of Perigee (rad)
    nu: float
        True Anomaly (rad)

    Notes
    -----
    This example is a real exercise from Orbital Mechanics for Engineering
    students by Howard D.Curtis. This exercise is 4.3 of 3rd. Edition, page 200.

    1. First the angular momentum is computed:

    .. math::
        \vec{h} = \vec{r} \times \vec{v}

    2. With it the eccentricity can be solved:

    .. math::
        \begin{align}
        \vec{e} &= \frac{1}{\mu}\left [ \left ( v^{2} - \frac{\mu}{r}\right ) \vec{r}  - (\vec{r} \cdot \vec{v})\vec{v} \right ] \\
        e &= \sqrt{\vec{e}\cdot\vec{e}} \\
        \end{align}

    3. The node vector line is solved:

    .. math::
        \begin{align}
        \vec{N} &= \vec{k} \times \vec{h} \\
        N &= \sqrt{\vec{N}\cdot\vec{N}}
        \end{align}

    4. The rigth ascension node is computed:

    .. math::
        \Omega = \left\{ \begin{array}{lcc}
         cos^{-1}{\left ( \frac{N_{x}}{N} \right )} &   if  & N_{y} \geq  0 \\
         \\ 360^{o} -cos^{-1}{\left ( \frac{N_{x}}{N} \right )} &  if & N_{y} < 0 \\
         \end{array}
        \right.

    5. The argument of perigee:

    .. math::
        \omega  = \left\{ \begin{array}{lcc}
         cos^{-1}{\left ( \frac{\vec{N}\vec{e}}{Ne} \right )} &   if  & e_{z} \geq  0 \\
         \\ 360^{o} -cos^{-1}{\left ( \frac{\vec{N}\vec{e}}{Ne} \right )} &  if & e_{z} < 0 \\
         \end{array}
        \right.

    6. And finally the true anomaly:

    .. math::
        \nu  = \left\{ \begin{array}{lcc}
         cos^{-1}{\left ( \frac{\vec{e}\vec{r}}{er} \right )} &   if  & v_{r} \geq  0 \\
         \\ 360^{o} -cos^{-1}{\left ( \frac{\vec{e}\vec{r}}{er} \right )} &  if & v_{r} < 0 \\
         \end{array}
        \right.

    Examples
    --------
    >>> from poliastro.bodies import Earth
    >>> from astropy import units as u
    >>> k = Earth.k.to_value(u.km ** 3 / u.s ** 2)
    >>> r = np.array([-6045., -3490., 2500.])
    >>> v = np.array([-3.457, 6.618, 2.533])
    >>> p, ecc, inc, raan, argp, nu = rv2coe(k, r, v)
    >>> print("p:", p, "[km]")  # doctest: +FLOAT_CMP
    p: 8530.47436396927 [km]
    >>> print("ecc:", ecc)  # doctest: +FLOAT_CMP
    ecc: 0.17121118195416898
    >>> print("inc:", np.rad2deg(inc), "[deg]")  # doctest: +FLOAT_CMP
    inc: 153.2492285182475 [deg]
    >>> print("raan:", np.rad2deg(raan), "[deg]")  # doctest: +FLOAT_CMP
    raan: 255.27928533439618 [deg]
    >>> print("argp:", np.rad2deg(argp), "[deg]")  # doctest: +FLOAT_CMP
    argp: 20.068139973005362 [deg]
    >>> print("nu:", np.rad2deg(nu), "[deg]")  # doctest: +FLOAT_CMP
    nu: 28.445804984192122 [deg]

    """

    h = cross(r, v)
    n = cross([0, 0, 1], h)
    e = ((v.dot(v) - k / (norm(r))) * r - r.dot(v) * v) / k
    ecc = norm(e)
    p = h.dot(h) / k
    inc = np.arccos(h[2] / norm(h))

    circular = ecc < tol
    equatorial = abs(inc) < tol

    if equatorial and not circular:
        raan = 0
        argp = np.arctan2(e[1], e[0]) % (2 * np.pi)  # Longitude of periapsis
        nu = np.arctan2(h.dot(cross(e, r)) / norm(h), r.dot(e))
    elif not equatorial and circular:
        raan = np.arctan2(n[1], n[0]) % (2 * np.pi)
        argp = 0
        # Argument of latitude
        nu = np.arctan2(r.dot(cross(h, n)) / norm(h), r.dot(n))
    elif equatorial and circular:
        raan = 0
        argp = 0
        nu = np.arctan2(r[1], r[0]) % (2 * np.pi)  # True longitude
    else:
        a = p / (1 - (ecc ** 2))
        ka = k * a
        if a > 0:
            e_se = r.dot(v) / sqrt(ka)
            e_ce = norm(r) * v.dot(v) / k - 1
            nu = E_to_nu(np.arctan2(e_se, e_ce), ecc)
        else:
            e_sh = r.dot(v) / sqrt(-ka)
            e_ch = norm(r) * (norm(v) ** 2) / k - 1
            nu = F_to_nu(np.log((e_ch + e_sh) / (e_ch - e_sh)) / 2, ecc)

        raan = np.arctan2(n[1], n[0]) % (2 * np.pi)
        px = r.dot(n)
        py = r.dot(cross(h, n)) / norm(h)
        argp = (np.arctan2(py, px) - nu) % (2 * np.pi)

    nu = (nu + np.pi) % (2 * np.pi) - np.pi

    return p, ecc, inc, raan, argp, nu


@jit
def mee2coe(p, f, g, h, k, L):
    r"""Converts from modified equinoctial orbital elements to classical
    orbital elements.

    The definition of the modified equinoctial orbital elements is taken from
    [Walker, 1985].

    .. math::

        \begin{align}
            p &= a(1 - e^{2})\\
            e &= \sqrt{f^{2} + g^{2}}\\
            i &= 2\arctan{(\sqrt{h^{2} + k^{2}})}\\
            raan &= atan2(k, h) \pmod{2\pi}\\
            argp &= (atan2(g, f) - raan) \pmod{2\pi}\\
            nu &= (L - atan2(g, f)) \pmod{2\pi}\\
        \end{align}

    Parameters
    ----------
    p: float
        Semi-latus rectum
    f: float
        Equinoctial parameter f
    g: float
        Equinoctial parameter g
    h: float
        Equinoctial parameter h
    k: float
        Equinoctial parameter k
    L: float
        Longitude

    Returns
    -------
    p: float
        Semi-latus rectum
    ecc: float
        Eccentricity of the orbit
    inc: float
        Inclination of the orbit
    raan: float
        RAAN of orbit
    argp: float
        Argument of the periapsis
    nu: float
        True anomaly

    Note
    -----
    The conversion is always safe because arctan2 works also for 0, 0
    arguments.

    """
    ecc = np.sqrt(f ** 2 + g ** 2)
    inc = 2 * np.arctan(np.sqrt(h ** 2 + k ** 2))
    lonper = np.arctan2(g, f)
    raan = np.arctan2(k, h) % (2 * np.pi)
    argp = (lonper - raan) % (2 * np.pi)
    nu = (L - lonper) % (2 * np.pi)
    return p, ecc, inc, raan, argp, nu
