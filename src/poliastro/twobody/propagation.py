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
    if ecc >= 1.0:
        raise ValueError("Mikkola only works with elliptical orbits.")
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


def _markley(k, r0, v0, tof):
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
    if ecc >= 1.0:
        raise ValueError("Markley only works with elliptical orbits.")
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


def pimienta(k, r, v, tofs, rtol=None):
    """ Kepler solver for both elliptic and parabolic orbits based on a 15th
    order polynomial with accuracies around 10e-5 for elliptic case and 10e-13
    in the hyperbolic regime.

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
    This algorithm was developed by Pimienta-Peñalver and John L. Crassidis in
    their paper *Accurate Kepler Equation solver without trascendental function
    evaluations*.
    """

    k = k.to(u.m ** 3 / u.s ** 2).value
    r0 = r.to(u.m).value
    v0 = v.to(u.m / u.s).value
    tofs = tofs.to(u.s).value

    results = [_pimienta(k, r0, v0, tof) for tof in tofs]
    return (
        [result[0] for result in results] * u.m,
        [result[1] for result in results] * u.m / u.s,
    )


def _pimienta(k, r0, v0, tof):
    """ Raw algorithm for Adonis' Pimienta and John L. Crassidis 15th order
    polynomial Kepler solver.

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
    """

    # Solve first for eccentricity and mean anomaly
    p, ecc, inc, raan, argp, nu = rv2coe(k, r0, v0)
    if ecc > 1.0:
        # TODO: include hyperbolic solver
        raise ValueError("Pimienta only works with elliptic orbits.")
    M0 = nu_to_M(nu, ecc)
    semi_axis_a = p / (1 - ecc ** 2)
    n = np.sqrt(k / np.abs(semi_axis_a) ** 3)
    M = M0 + n * tof

    # Equation (32a), (32b), (32c) and (32d)
    c3 = 5 / 2 + 560 * ecc
    a = 15 * (1 - ecc) / c3
    b = -M / c3
    y = np.sqrt(b ** 2 / 4 + a ** 3 / 27)

    # Equation (33)
    x_bar = (-b / 2 + y) ** (1 / 3) - (b / 2 + y) ** (1 / 3)

    # Coefficients from equations (34a) and (34b)
    c15 = 3003 / 14336 + 16384 * ecc
    c13 = 3465 / 13312 - 61440 * ecc
    c11 = 945 / 2816 + 92160 * ecc
    c9 = 175 / 384 - 70400 * ecc
    c7 = 75 / 112 + 28800 * ecc
    c5 = 9 / 8 - 6048 * ecc

    # Precompute x_bar powers, equations (35a) to (35d)
    x_bar2 = x_bar ** 2
    x_bar3 = x_bar2 * x_bar
    x_bar4 = x_bar3 * x_bar
    x_bar5 = x_bar4 * x_bar
    x_bar6 = x_bar5 * x_bar
    x_bar7 = x_bar6 * x_bar
    x_bar8 = x_bar7 * x_bar
    x_bar9 = x_bar8 * x_bar
    x_bar10 = x_bar9 * x_bar
    x_bar11 = x_bar10 * x_bar
    x_bar12 = x_bar11 * x_bar
    x_bar13 = x_bar12 * x_bar
    x_bar14 = x_bar13 * x_bar
    x_bar15 = x_bar14 * x_bar

    # Function f and its derivatives are given by all the (36) equation set
    f = (
        c15 * x_bar15
        + c13 * x_bar13
        + c11 * x_bar11
        + c9 * x_bar9
        + c7 * x_bar7
        + c5 * x_bar5
        + c3 * x_bar3
        + 15 * (1 - ecc) * x_bar
        - M
    )
    f1 = (
        15 * c15 * x_bar14
        + 13 * c13 * x_bar12
        + 11 * c11 * x_bar10
        + 9 * c9 * x_bar8
        + 7 * c7 * x_bar6
        + 5 * c5 * x_bar4
        + 3 * c3 * x_bar2
        + 15 * (1 - ecc)
    )
    f2 = (
        210 * c15 * x_bar13
        + 156 * c13 * x_bar11
        + 110 * c11 * x_bar9
        + 72 * c9 * x_bar7
        + 42 * c7 * x_bar5
        + 20 * c5 * x_bar3
        + 6 * c3 * x_bar
    )
    f3 = (
        2730 * c15 * x_bar12
        + 1716 * c13 * x_bar10
        + 990 * c11 * x_bar8
        + 504 * c9 * x_bar6
        + 210 * c7 * x_bar4
        + 60 * c5 * x_bar2
        + 6 * c3
    )
    f4 = (
        32760 * c15 * x_bar11
        + 17160 * c13 * x_bar9
        + 7920 * c11 * x_bar7
        + 3024 * c9 * x_bar5
        + 840 * c7 * x_bar3
        + 120 * c5 * x_bar
    )
    f5 = (
        360360 * c15 * x_bar10
        + 154440 * c13 * x_bar8
        + 55440 * c11 * x_bar6
        + 15120 * c9 * x_bar4
        + 2520 * c7 * x_bar2
        + 120 * c5
    )
    f6 = (
        3603600 * c15 * x_bar9
        + 1235520 * c13 * x_bar7
        + 332640 * c11 * x_bar5
        + 60480 * c9 * x_bar3
        + 5040 * c7 * x_bar
    )
    f7 = (
        32432400 * c15 * x_bar8
        + 8648640 * c13 * x_bar6
        + 1663200 * c11 * x_bar4
        + 181440 * c9 * x_bar2
        + 5040 * c7
    )
    f8 = (
        259459200 * c15 * x_bar7
        + 51891840 * c13 * x_bar5
        + 6652800 * c11 * x_bar3
        + 362880 * c9 * x_bar
    )
    f9 = (
        1.8162144e9 * c15 * x_bar6
        + 259459200 * c13 * x_bar4
        + 19958400 * c11 * x_bar2
        + 362880 * c9
    )
    f10 = (
        1.08972864e10 * c15 * x_bar5
        + 1.0378368e9 * c13 * x_bar3
        + 39916800 * c11 * x_bar
    )
    f11 = 5.4486432e10 * c15 * x_bar4 + 3.1135104e9 * c13 * x_bar2 + 39916800 * c11
    f12 = 2.17945728e11 * c15 * x_bar3 + 6.2270208e9 * c13 * x_bar
    f13 = 6.53837184 * c15 * x_bar2 + 6.2270208e9 * c13
    f14 = 1.307674368e13 * c15 * x_bar
    f15 = 1.307674368e13 * c15

    # Solving g parameters defined by equations (37a), (37b), (37c) and (37d)
    g1 = 1 / 2
    g2 = 1 / 6
    g3 = 1 / 24
    g4 = 1 / 120
    g5 = 1 / 720
    g6 = 1 / 5040
    g7 = 1 / 40320
    g8 = 1 / 362880
    g9 = 1 / 3628800
    g10 = 1 / 39916800
    g11 = 1 / 479001600
    g12 = 1 / 6.2270208e9
    g13 = 1 / 8.71782912e10
    g14 = 1 / 1.307674368e12

    # Solving for the u_{i} and h_{i} variables defined by equation (38)
    u1 = -f / f1

    h2 = f1 + g1 * u1 * f2
    u2 = -f / h2

    h3 = f1 + g1 * u2 * f2 + g2 * u2 ** 2 * f3
    u3 = -f / h3

    h4 = f1 + g1 * u3 * f2 + g2 * u3 ** 2 * f3 + g3 * u3 ** 3 * f4
    u4 = -f / h4

    h5 = f1 + g1 * u4 * f2 + g2 * u4 ** 2 * f3 + g3 * u4 ** 3 * f4 + g4 * u4 ** 4 * f5
    u5 = -f / h5

    h6 = (
        f1
        + g1 * u5 * f2
        + g2 * u5 ** 2 * f3
        + g3 * u5 ** 3 * f4
        + g4 * u5 ** 4 * f5
        + g5 * u5 ** 5 * f6
    )
    u6 = -f / h6

    h7 = (
        f1
        + g1 * u6 * f2
        + g2 * u6 ** 2 * f3
        + g3 * u6 ** 3 * f4
        + g4 * u6 ** 4 * f5
        + g5 * u6 ** 5 * f6
        + g6 * u6 ** 6 * f7
    )
    u7 = -f / h7

    h8 = (
        f1
        + g1 * u7 * f2
        + g2 * u7 ** 2 * f3
        + g3 * u7 ** 3 * f4
        + g4 * u7 ** 4 * f5
        + g5 * u7 ** 5 * f6
        + g6 * u7 ** 6 * f7
        + g7 * u7 ** 7 * f8
    )
    u8 = -f / h8

    h9 = (
        f1
        + g1 * u8 * f2
        + g2 * u8 ** 2 * f3
        + g3 * u8 ** 3 * f4
        + g4 * u8 ** 4 * f5
        + g5 * u8 ** 5 * f6
        + g6 * u8 ** 6 * f7
        + g7 * u8 ** 7 * f8
        + g8 * u8 ** 8 * f9
    )
    u9 = -f / h9

    h10 = (
        f1
        + g1 * u9 * f2
        + g2 * u9 ** 2 * f3
        + g3 * u9 ** 3 * f4
        + g4 * u9 ** 4 * f5
        + g5 * u9 ** 5 * f6
        + g6 * u9 ** 6 * f7
        + g7 * u9 ** 7 * f8
        + g8 * u9 ** 8 * f9
        + g9 * u9 ** 9 * f10
    )
    u10 = -f / h10

    h11 = (
        f1
        + g1 * u10 * f2
        + g2 * u10 ** 2 * f3
        + g3 * u10 ** 3 * f4
        + g4 * u10 ** 4 * f5
        + g5 * u10 ** 5 * f6
        + g6 * u10 ** 6 * f7
        + g7 * u10 ** 7 * f8
        + g8 * u10 ** 8 * f9
        + g9 * u10 ** 9 * f10
        + g10 * u10 ** 10 * f11
    )
    u11 = -f / h11

    h12 = (
        f1
        + g1 * u11 * f2
        + g2 * u11 ** 2 * f3
        + g3 * u11 ** 3 * f4
        + g4 * u11 ** 4 * f5
        + g5 * u11 ** 5 * f6
        + g6 * u11 ** 6 * f7
        + g7 * u11 ** 7 * f8
        + g8 * u11 ** 8 * f9
        + g9 * u11 ** 9 * f10
        + g10 * u11 ** 10 * f11
        + g11 * u11 ** 11 * f12
    )
    u12 = -f / h12

    h13 = (
        f1
        + g1 * u12 * f2
        + g2 * u12 ** 2 * f3
        + g3 * u12 ** 3 * f4
        + g4 * u12 ** 4 * f5
        + g5 * u12 ** 5 * f6
        + g6 * u12 ** 6 * f7
        + g7 * u12 ** 7 * f8
        + g8 * u12 ** 8 * f9
        + g9 * u12 ** 9 * f10
        + g10 * u12 ** 10 * f11
        + g11 * u12 ** 11 * f12
        + g12 * u12 ** 12 * f13
    )
    u13 = -f / h13

    h14 = (
        f1
        + g1 * u13 * f2
        + g2 * u13 ** 2 * f3
        + g3 * u13 ** 3 * f4
        + g4 * u13 ** 4 * f5
        + g5 * u13 ** 5 * f6
        + g6 * u13 ** 6 * f7
        + g7 * u13 ** 7 * f8
        + g8 * u13 ** 8 * f9
        + g9 * u13 ** 9 * f10
        + g10 * u13 ** 10 * f11
        + g11 * u13 ** 11 * f12
        + g12 * u13 ** 12 * f13
        + g13 * u13 ** 13 * f14
    )
    u14 = -f / h14

    h15 = (
        f1
        + g1 * u14 * f2
        + g2 * u14 ** 2 * f3
        + g3 * u14 ** 3 * f4
        + g4 * u14 ** 4 * f5
        + g5 * u14 ** 5 * f6
        + g6 * u14 ** 6 * f7
        + g7 * u14 ** 7 * f8
        + g8 * u14 ** 8 * f9
        + g9 * u14 ** 9 * f10
        + g10 * u14 ** 10 * f11
        + g11 * u14 ** 11 * f12
        + g12 * u14 ** 12 * f13
        + g13 * u14 ** 13 * f14
        + g14 * u14 ** 14 * f15
    )
    u15 = -f / h15

    # Solving for x
    x = x_bar + u15
    w = x - 0.01171875 * x ** 17 / (1 + ecc)

    # Solving for the true anomaly from eccentricity anomaly
    E = M + ecc * (
        -16384 * w ** 15
        + 61440 * w ** 13
        - 92160 * w ** 11
        + 70400 * w ** 9
        - 28800 * w ** 7
        + 6048 * w ** 5
        - 560 * w ** 3
        + 15 * w
    )

    nu = E_to_nu(E, ecc)

    return coe2rv(k, p, ecc, inc, raan, argp, nu)


def goodingI(k, r, v, tofs, numiter=150, rtol=1e-8):
    """ Solves the Elliptic Kepler Equation with a cubic convergence and
    accuracy better than 10e-12 rad is normally achieved. It is not accurate
    when ecc=1.0 and M=0.

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
    This method was developed by Gooding and Odell in their paper *The
    hyperbolic Kepler equation (and the elliptic equation revisited)* with
    DOI: https://doi.org/10.1007/BF01235540
    """

    k = k.to(u.m ** 3 / u.s ** 2).value
    r0 = r.to(u.m).value
    v0 = v.to(u.m / u.s).value
    tofs = tofs.to(u.s).value

    results = [_goodingI(k, r0, v0, tof, numiter=numiter, rtol=rtol) for tof in tofs]
    return (
        [result[0] for result in results] * u.m,
        [result[1] for result in results] * u.m / u.s,
    )


def _goodingI(k, r0, v0, tof, numiter=150, rtol=1e-8):
    """ Solves the Elliptic Kepler Equation with a cubic convergence and
    accuracy better than 10e-12 rad is normally achieved. It is not accurate
    when ecc=1.0 and M=0.

    Parameters
    ----------
    k : float
        Standard gravitational parameter of the attractor.
    r : 1x3 vector
        Position vector.
    v : 1x3 vector
        Velocity vector.
    tof : float
        Time of flight.
    rtol: float
        Relative error for accuracy of the method.

    Returns
    -------
    rr : 1x3 vector
        Propagated position vectors.
     vv : 1x3 vector

    Note
    ----
    Original paper for the algorithm: https://doi.org/10.1007/BF01238923
    """

    # Solve first for eccentricity and mean anomaly
    p, ecc, inc, raan, argp, nu = rv2coe(k, r0, v0)
    if ecc >= 1.0:
        raise ValueError("Gooding I only works with elliptic orbits.")
    M0 = nu_to_M(nu, ecc)
    semi_axis_a = p / (1 - ecc ** 2)
    n = np.sqrt(k / np.abs(semi_axis_a) ** 3)
    M = M0 + n * tof

    # Start the computation
    n = 0
    c = ecc * np.cos(M)
    s = ecc * np.sin(M)
    psi = s / np.sqrt(1 - 2 * c + ecc ** 2)
    f = 1.0
    while f ** 2 >= rtol and n <= numiter:
        xi = np.cos(psi)
        eta = np.sin(psi)
        fd = (1 - c * xi) + s * eta
        fdd = c * eta + s * xi
        f = psi - fdd
        psi = psi - f * fd / (fd ** 2 - 0.5 * f * fdd)
        n += 1

    E = M + psi
    nu = E_to_nu(E, ecc)

    return coe2rv(k, p, ecc, inc, raan, argp, nu)


def danby(k, r, v, tofs, rtol=1e-8):
    """ Kepler solver for both elliptic and parabolic orbits based on Danby's
    algorithm.

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
        Relative error for accuracy of the method.

    Returns
    -------
    rr : ~astropy.units.Quantity
        Propagated position vectors.
    vv : ~astropy.units.Quantity
        Propagated velocity vectors.

    Note
    ----
    This algorithm was developed by Danby in his paper *The solution of Kepler
    Equation* with DOI: https://doi.org/10.1007/BF01686811
    """

    k = k.to(u.m ** 3 / u.s ** 2).value
    r0 = r.to(u.m).value
    v0 = v.to(u.m / u.s).value
    tofs = tofs.to(u.s).value

    results = [_danby(k, r0, v0, tof) for tof in tofs]
    return (
        [result[0] for result in results] * u.m,
        [result[1] for result in results] * u.m / u.s,
    )


def _danby(k, r0, v0, tof, numiter=20, rtol=1e-8):
    """ Kepler solver for both elliptic and parabolic orbits based on Danby's
    algorithm.

    Parameters
    ----------
    k : float
        Standard gravitational parameter of the attractor.
    r : 1x3 vector
        Position vector.
    v : 1x3 vector
        Velocity vector.
    tof : float
        Time of flight.
    rtol: float
        Relative error for accuracy of the method.

    Returns
    -------
    rr : 1x3 vector
        Propagated position vectors.
    vv : 1x3 vector

    Note
    ----
    This algorithm was developed by Danby in his paper *The solution of Kepler
    Equation* with DOI: https://doi.org/10.1007/BF01686811
    """

    # Solve first for eccentricity and mean anomaly
    p, ecc, inc, raan, argp, nu = rv2coe(k, r0, v0)
    M0 = nu_to_M(nu, ecc)
    semi_axis_a = p / (1 - ecc ** 2)
    n = np.sqrt(k / np.abs(semi_axis_a) ** 3)
    M = M0 + n * tof

    # Range mean anomaly
    xma = M - 2 * np.pi * np.floor(M / 2 / np.pi)

    if ecc == 0:
        # Solving for circular orbit
        nu = xma
        return coe2rv(k, p, ecc, inc, raan, argp, nu)

    elif ecc < 1.0:
        # For elliptical orbit
        E = xma + 0.85 * np.sign(np.sin(xma)) * ecc

    else:
        # For parabolic and hyperbolic
        E = np.log(2 * xma / ecc + 1.8)

    # Iterations begin
    n = 0
    while n <= numiter:

        if ecc < 1.0:
            s = ecc * np.sin(E)
            c = ecc * np.cos(E)
            f = E - s - xma
            fp = 1 - c
            fpp = s
            fppp = c
        else:
            s = ecc * np.sinh(E)
            c = ecc * np.cosh(E)
            f = s - E - xma
            fp = c - 1
            fpp = s
            fppp = c

        if np.abs(f) <= rtol:

            if ecc < 1.0:
                sta = np.sqrt(1 - ecc ** 2) * np.sin(E)
                cta = np.cos(E) - ecc
            else:
                sta = np.sqrt(ecc ** 2 - 1) * np.sinh(E)
                cta = ecc - np.cosh(E)

            nu = np.arctan2(sta, cta)
            return coe2rv(k, p, ecc, inc, raan, argp, nu)

        else:
            delta = -f / fp
            delta_star = -f / (fp + 0.5 * delta * fpp)
            deltak = -f / (fp + 0.5 * delta_star * fpp + delta_star ** 2 * fppp / 6)
            E = E + deltak
            n += 1

    raise ValueError("Maximum number of iterations has been reached.")


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
