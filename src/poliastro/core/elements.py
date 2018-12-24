"""This module contains a set of functions that can be used
    to convert between different elements that define the orbit
    of a body.
    """

import numpy as np
from numpy.core.umath import cos, sin, sqrt

from poliastro.core.util import cross, norm, transform

from ._jit import jit


@jit
def rv_pqw(k, p, ecc, nu):
    r"""Returns r and v vectors in perifocal frame.

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

    Examples
    --------
    >>> from poliastro.core.elements import rv_pqw
    >>> from poliastro.constants import GM_earth
    >>> k = GM_earth #Earth gravitational parameter

    >>> ecc = 0.3 #Eccentricity
    >>> h = 60000e6 #Angular momentum of the orbit [m^2]/[s]
    >>> nu = np.deg2rad(120) #True Anomaly [rad]
    >>> p = h**2 / k #Parameter of the orbit
    >>> r, v = rv_pqw(k, p, ecc, nu)

    >>> #Printing the results
    r = [-5312706.25105345  9201877.15251336    0] [m]
    v = [-5753.30180931 -1328.66813933  0] [m]/[s]

    Note
    ----

    These formulas can be checked at Curtis 3rd. Edition, page 110. Also the example proposed is 2.11
    of Curtis 3rd Edition book.
    """
    r_pqw = (np.array([cos(nu), sin(nu), 0 * nu]) * p / (1 + ecc * cos(nu))).T
    v_pqw = (np.array([-sin(nu), (ecc + cos(nu)), 0]) * sqrt(k / p)).T
    return r_pqw, v_pqw


@jit
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

    r_ijk = transform(r_pqw, -argp, 2)
    r_ijk = transform(r_ijk, -inc, 0)
    r_ijk = transform(r_ijk, -raan, 2)
    v_ijk = transform(v_pqw, -argp, 2)
    v_ijk = transform(v_ijk, -inc, 0)
    v_ijk = transform(v_ijk, -raan, 2)

    return r_ijk, v_ijk


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

    .. math::
        \begin{align}
        p &= a(1-e^2) \\
        f &= e\cos(\omega + \Omega) \\
        g &= e\sin(\omega + \Omega) \\
        h &= \tan(\frac{i}{2})\cos(\Omega) \\
        k &= \tan(\frac{i}{2})\sin(\Omega) \\
        L &= \Omega + \omega + \theta \\
        \end{align}
    
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

    Note
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


@jit
def rv2coe(k, r, v, tol=1e-8):
    r"""Converts from vectors to classical orbital elements.

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

    Examples
    --------
    >>> from poliastro.core.elements import rv2coe
    >>> from poliastro.constants import GM_earth
    >>> from astropy import units as u
    >>> import numpy as np
    >>> k = GM_earth.to(u.km**3 / u.s**2) #Earth gravitational parameter
    >>> r = [-6045, -3490, 2500] * u.km
    >>> v = [-3.457, 6.618, 2.533] * u.km / u.s
    >>> p, ecc, inc, raan, argp, nu = rv2coe(k, r, v)
    >>> print("p:", p, "[km]")
    p: 8530.47436396927 [km]
    >>> print("ecc:", ecc)
    ecc: 0.17121118195416898
    >>> print("inc:", np.rad2deg(inc), "[deg]")
    inc: 153.2492285182475 [deg]
    >>> print("raan:", np.rad2deg(raan), "[deg]")
    raan: 255.27928533439618 [deg]
    >>> print("argp:", np.rad2deg(argp), "[deg]")
    argp: 20.068139973005362 [deg]
    >>> print("nu:", np.rad2deg(nu), "[deg]")
    nu: 28.445804984192108 [deg]


    Note
    ----

    This example is a real exercise from Orbital Mechanics for Engineering
    students by Howard D.Curtis. This exercise is 4.3 of 3rd. Edition, page 200.
    """

    h = cross(r, v)
    n = cross([0, 0, 1], h) / norm(h)
    e = ((v.dot(v) - k / (norm(r))) * r - r.dot(v) * v) / k
    ecc = norm(e)
    p = h.dot(h) / k
    inc = np.arccos(h[2] / norm(h))

    circular = ecc < tol
    equatorial = abs(inc) < tol

    if equatorial and not circular:
        raan = 0
        argp = np.arctan2(e[1], e[0]) % (2 * np.pi)  # Longitude of periapsis
        nu = np.arctan2(h.dot(cross(e, r)) / norm(h), r.dot(e)) % (2 * np.pi)
    elif not equatorial and circular:
        raan = np.arctan2(n[1], n[0]) % (2 * np.pi)
        argp = 0
        # Argument of latitude
        nu = np.arctan2(r.dot(cross(h, n)) / norm(h), r.dot(n)) % (2 * np.pi)
    elif equatorial and circular:
        raan = 0
        argp = 0
        nu = np.arctan2(r[1], r[0]) % (2 * np.pi)  # True longitude
    else:
        raan = np.arctan2(n[1], n[0]) % (2 * np.pi)
        argp = np.arctan2(e.dot(cross(h, n)) / norm(h), e.dot(n)) % (2 * np.pi)
        nu = np.arctan2(r.dot(cross(h, e)) / norm(h), r.dot(e)) % (2 * np.pi)

    return p, ecc, inc, raan, argp, nu
