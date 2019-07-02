"""Propagation algorithms

"""
import functools

import numpy as np
from astropy import time, units as u
from astropy.coordinates import CartesianDifferential, CartesianRepresentation
from scipy.integrate import solve_ivp

from poliastro.core.angles import (
    E_to_nu,
    _kepler_equation,
    _kepler_equation_prime,
    nu_to_M,
    F_to_nu
)
from poliastro.core.elements import coe2rv, rv2coe
from poliastro.core.propagation import (
    func_twobody,
    kepler as kepler_fast,
    mean_motion as mean_motion_fast,
)
from poliastro.integrators import DOP835


def cowell(k, r, v, tofs, rtol=1e-11, *, ad=None, **ad_kwargs):
    """Propagates orbit using Cowell's formulation.

    Parameters
    ----------
    k : ~astropy.units.Quantity
        Standard gravitational parameter of the attractor.
    r : ~astropy.units.Quantity
        Position vector.
    v : ~astropy.units.Quantity
        Velocity vector.
    tofs : ~astropy.units.Quantity
        Array of times to propagate.
    rtol : float, optional
        Maximum relative error permitted, default to 1e-10.
    ad : function(t0, u, k), optional
         Non Keplerian acceleration (km/s2), default to None.

    Returns
    -------
    rr : ~astropy.units.Quantity
        Propagated position vectors.
    vv : ~astropy.units.Quantity
        Propagated velocity vectors.

    Raises
    ------
    RuntimeError
        If the algorithm didn't converge.

    Note
    -----
    This method uses a Dormand & Prince method of order 8(5,3) available
    in the :py:class:`poliastro.integrators` module. If multiple tofs
    are provided, the method propagates to the maximum value and
    calculates the other values via dense output

    """
    k = k.to(u.km ** 3 / u.s ** 2).value
    x, y, z = r.to(u.km).value
    vx, vy, vz = v.to(u.km / u.s).value
    tofs = tofs.to(u.s).value

    u0 = np.array([x, y, z, vx, vy, vz])

    # Set the non Keplerian acceleration
    if ad is None:

        def ad(t0, u_, k_):
            return 0, 0, 0

    f_with_ad = functools.partial(func_twobody, k=k, ad=ad, ad_kwargs=ad_kwargs)

    result = solve_ivp(
        f_with_ad,
        (0, max(tofs)),
        u0,
        rtol=rtol,
        atol=1e-12,
        method=DOP835,
        dense_output=True,
    )
    if not result.success:
        raise RuntimeError("Integration failed")

    rrs = []
    vvs = []
    for i in range(len(tofs)):
        y = result.sol(tofs[i])
        rrs.append(y[:3])
        vvs.append(y[3:])

    return rrs * u.km, vvs * u.km / u.s


def mean_motion(k, r, v, tofs, **kwargs):
    """Propagates orbit using Cowell's formulation.

    Parameters
    ----------
    k : ~astropy.units.Quantity
        Standard gravitational parameter of the attractor.
    r : ~astropy.units.Quantity
        Position vector.
    v : ~astropy.units.Quantity
        Velocity vector.
    tofs : ~astropy.units.Quantity
        Array of times to propagate.

    Returns
    -------
    rr : ~astropy.units.Quantity
        Propagated position vectors.
    vv : ~astropy.units.Quantity
        Propagated velocity vectors.

    """
    k = k.to(u.km ** 3 / u.s ** 2).value
    r0 = r.to(u.km).value
    v0 = v.to(u.km / u.s).value
    tofs = tofs.to(u.s).value

    results = [mean_motion_fast(k, r0, v0, tof) for tof in tofs]
    # TODO: Rewrite to avoid iterating twice
    return (
        [result[0] for result in results] * u.km,
        [result[1] for result in results] * u.km / u.s,
    )


def kepler(k, r, v, tofs, numiter=350, **kwargs):
    """Propagates Keplerian orbit.

    Parameters
    ----------
    k : ~astropy.units.Quantity
        Standard gravitational parameter of the attractor.
    r : ~astropy.units.Quantity
        Position vector.
    v : ~astropy.units.Quantity
        Velocity vector.
    tofs : ~astropy.units.Quantity
        Array of times to propagate.
    numiter : int, optional
        Maximum number of iterations, default to 35.

    Returns
    -------
    rr : ~astropy.units.Quantity
        Propagated position vectors.
    vv : ~astropy.units.Quantity
        Propagated velocity vectors.

    Raises
    ------
    RuntimeError
        If the algorithm didn't converge.

    Note
    -----
    This algorithm is based on Vallado implementation, and does basic Newton
    iteration on the Kepler equation written using universal variables. Battin
    claims his algorithm uses the same amount of memory but is between 40 %
    and 85 % faster.

    """
    k = k.to(u.km ** 3 / u.s ** 2).value
    r0 = r.to(u.km).value
    v0 = v.to(u.km / u.s).value
    tofs = tofs.to(u.s).value

    results = [_kepler(k, r0, v0, tof, numiter=numiter) for tof in tofs]
    # TODO: Rewrite to avoid iterating twice
    return (
        [result[0] for result in results] * u.km,
        [result[1] for result in results] * u.km / u.s,
    )


def _kepler(k, r0, v0, tof, *, numiter):
    # Compute Lagrange coefficients
    f, g, fdot, gdot = kepler_fast(k, r0, v0, tof, numiter)

    assert np.abs(f * gdot - fdot * g - 1) < 1e-5  # Fixed tolerance

    # Return position and velocity vectors
    r = f * r0 + g * v0
    v = fdot * r0 + gdot * v0

    return r, v


def mikkola(k, r, v, tofs, rtol=None):
    """ Solves Kepler Equation by a cubic approximation. This method is valid
    no mater the orbit's nature.

    Parameters
    ----------
    k : ~astropy.units.Quantity
        Standard gravitational parameter of the attractor.
    r : ~astropy.units.Quantity
        Position vector.
    v : ~astropy.units.Quantity
        Velocity vector.
    tofs : ~astropy.units.Quantity
        Array of times to propagate.
    rtol: float
        This method does not require of tolerance since it is non iterative.

    Returns
    -------
    rr : ~astropy.units.Quantity
        Propagated position vectors.
    vv : ~astropy.units.Quantity

    Note
    ----
    This method was derived by Seppo Mikola in his paper *A Cubic Approximation
    For Kepler's Equation* with DOI: https://doi.org/10.1007/BF01235850
    """

    k = k.to(u.m ** 3 / u.s ** 2).value
    r0 = r.to(u.m).value
    v0 = v.to(u.m / u.s).value
    tofs = tofs.to(u.s).value

    results = [_mikkola(k, r0, v0, tof) for tof in tofs]
    return (
        [result[0] for result in results] * u.m,
        [result[1] for result in results] * u.m / u.s,
    )


def _mikkola(k, r0, v0, tof, rtol=None):
    """ Raw algorithm for Mikkola's Kepler solver.

    Parameters
    ----------
    k : ~astropy.units.Quantity
        Standard gravitational parameter of the attractor.
    r : ~astropy.units.Quantity
        Position vector.
    v : ~astropy.units.Quantity
        Velocity vector.
    tofs : ~astropy.units.Quantity
        Array of times to propagate.
    rtol: float
        This method does not require of tolerance since it is non iterative.

    Returns
    -------
    rr : ~astropy.units.Quantity
        Propagated position vectors.
    vv : ~astropy.units.Quantity

    Note
    ----
    Original paper: https://doi.org/10.1007/BF01235850
    """

    # Solving for the classical elements
    p, ecc, inc, raan, argp, nu = rv2coe(k, r0, v0)
    M0 = nu_to_M(nu, ecc)
    a = p / (1 - ecc ** 2)
    n = np.sqrt(k / a ** 3)
    M = M0 + n * tof

    if ecc < 1.0:
        # Solving elliptical case

        # Equation (9a)
        alpha = (1 - ecc) / (4 * ecc + 1 / 2)
        beta = M / 2 / (4 * ecc + 1 / 2)

        # Equation (9b)
        if beta >= 0:
            z = (beta + np.sqrt(beta ** 2 + alpha ** 3)) ** (1 / 3)
        else:
            z = (beta - np.sqrt(beta ** 2 + alpha ** 3)) ** (1 / 3)

        # Equation (9c)
        s = z - alpha / z

        # Apply correction term
        ds = -0.078 * s ** 5 / (1 + ecc)
        s = s + ds

        # Equation (8)
        E = M + ecc * (3 * s - 4 * s ** 3)
        nu = E_to_nu(E, ecc)

        return coe2rv(k, p, ecc, inc, raan, argp, nu)

def markley(k, r, v, tofs, rtol=None):
    """ Kepler Equation solver based on a fifth-order refinement of the solution
    of a cubic equation.

    Parameters
    ----------
    k : ~astropy.units.Quantity
        Standard gravitational parameter of the attractor.
    r : ~astropy.units.Quantity
        Position vector.
    v : ~astropy.units.Quantity
        Velocity vector.
    tofs : ~astropy.units.Quantity
        Array of times to propagate.
    rtol: float
        This method does not require of tolerance since it is non iterative.

    Returns
    -------
    rr : ~astropy.units.Quantity
        Propagated position vectors.
    vv : ~astropy.units.Quantity
        Propagated velocity vectors.

    Note
    ----
    This method was originally presented by Markley in his paper *Kepler Equation Solver*
    with DOI: https://doi.org/10.1007/BF00691917
    """

    k = k.to(u.m ** 3 / u.s ** 2).value
    r0 = r.to(u.m).value
    v0 = v.to(u.m / u.s).value
    tofs = tofs.to(u.s).value

    results = [_markley(k, r0, v0, tof) for tof in tofs]
    return (
        [result[0] for result in results] * u.m,
        [result[1] for result in results] * u.m / u.s,
    )


def _markley(k, r0, v0, tof, rtol=None):
    """ Solves the kepler problem by a non iterative method. Relative error is
    around 1e-18, only limited by machine double-precission errors.

    Parameters
    ----------
    k : float
        Standar Gravitational parameter
    r0 : array
        Initial position vector wrt attractor center.
    v0 : array
        Initial velocity vector.
    tof : float
        Time of flight.

    Returns
    -------
    rf: array
        Final position vector
    vf: array
        Final velocity vector

    Note
    ----
    The following algorithm was taken from http://dx.doi.org/10.1007/BF00691917.
    """

    # Solve first for eccentricity and mean anomaly
    p, ecc, inc, raan, argp, nu = rv2coe(k, r0, v0)
    M0 = nu_to_M(nu, ecc)
    a = p / (1 - ecc ** 2)
    n = np.sqrt(k / a ** 3)
    M = M0 + n * tof

    # Range between -pi and pi
    M = M % (2 * np.pi)
    if M > np.pi:
        M = -(2 * np.pi - M)

    # Equation (20)
    alpha = (3 * np.pi ** 2 + 1.6 * (np.pi - np.abs(M)) / (1 + ecc)) / (np.pi ** 2 - 6)

    # Equation (5)
    d = 3 * (1 - ecc) + alpha * ecc

    # Equation (9)
    q = 2 * alpha * d * (1 - ecc) - M ** 2

    # Equation (10)
    r = 3 * alpha * d * (d - 1 + ecc) * M + M ** 3

    # Equation (14)
    w = (np.abs(r) + np.sqrt(q ** 3 + r ** 2)) ** (2 / 3)

    # Equation (15)
    E = (2 * r * w / (w ** 2 + w * q + q ** 2) + M) / d

    # Equation (26)
    f0 = _kepler_equation(E, M, ecc)
    f1 = _kepler_equation_prime(E, M, ecc)
    f2 = ecc * np.sin(E)
    f3 = ecc * np.cos(E)
    f4 = -f2

    # Equation (22)
    delta3 = -f0 / (f1 - 0.5 * f0 * f2 / f1)
    delta4 = -f0 / (f1 + 0.5 * delta3 * f2 + 1 / 6 * delta3 ** 2 * f3)
    delta5 = -f0 / (
        f1 + 0.5 * delta4 * f2 + 1 / 6 * delta4 ** 2 * f3 + 1 / 24 * delta4 ** 3 * f4
    )

    E += delta5
    nu = E_to_nu(E, ecc)

    return coe2rv(k, p, ecc, inc, raan, argp, nu)


def propagate(orbit, time_of_flight, *, method=mean_motion, rtol=1e-10, **kwargs):
    """Propagate an orbit some time and return the result.

    Parameters
    ----------
    orbit : ~poliastro.twobody.Orbit
        Orbit object to propagate.
    time_of_flight : ~astropy.time.TimeDelta
        Time of propagation.
    method : callable, optional
        Propagation method, default to mean_motion.
    rtol : float, optional
        Relative tolerance, default to 1e-10.

    """
    # Use the highest precision we can afford
    # np.atleast_1d does not work directly on TimeDelta objects
    jd1 = np.atleast_1d(time_of_flight.jd1)
    jd2 = np.atleast_1d(time_of_flight.jd2)
    time_of_flight = time.TimeDelta(jd1, jd2, format="jd", scale=time_of_flight.scale)

    rr, vv = method(
        orbit.attractor.k, orbit.r, orbit.v, time_of_flight.to(u.s), rtol=rtol, **kwargs
    )

    # TODO: Turn these into unit tests
    assert rr.ndim == 2
    assert vv.ndim == 2

    cartesian = CartesianRepresentation(
        rr, differentials=CartesianDifferential(vv, xyz_axis=1), xyz_axis=1
    )

    return cartesian
