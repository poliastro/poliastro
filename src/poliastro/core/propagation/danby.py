import numpy as np
from numba import njit as jit

from poliastro.core.angles import E_to_M, F_to_M, nu_to_E, nu_to_F
from poliastro.core.elements import coe2rv, rv2coe


@jit
def danby_coe(k, p, ecc, inc, raan, argp, nu, tof, numiter=20, rtol=1e-8):

    semi_axis_a = p / (1 - ecc**2)
    n = np.sqrt(k / np.abs(semi_axis_a) ** 3)

    if ecc == 0:
        # Solving for circular orbit
        M0 = nu  # for circular orbit M = E = nu
        M = M0 + n * tof
        nu = M - 2 * np.pi * np.floor(M / 2 / np.pi)
        return nu

    elif ecc < 1.0:
        # For elliptical orbit
        M0 = E_to_M(nu_to_E(nu, ecc), ecc)
        M = M0 + n * tof
        xma = M - 2 * np.pi * np.floor(M / 2 / np.pi)
        E = xma + 0.85 * np.sign(np.sin(xma)) * ecc

    else:
        # For parabolic and hyperbolic
        M0 = F_to_M(nu_to_F(nu, ecc), ecc)
        M = M0 + n * tof
        xma = M - 2 * np.pi * np.floor(M / 2 / np.pi)
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
                sta = np.sqrt(1 - ecc**2) * np.sin(E)
                cta = np.cos(E) - ecc
            else:
                sta = np.sqrt(ecc**2 - 1) * np.sinh(E)
                cta = ecc - np.cosh(E)

            nu = np.arctan2(sta, cta)
            break
        else:
            delta = -f / fp
            delta_star = -f / (fp + 0.5 * delta * fpp)
            deltak = -f / (
                fp + 0.5 * delta_star * fpp + delta_star**2 * fppp / 6
            )
            E = E + deltak
            n += 1
    else:
        raise ValueError("Maximum number of iterations has been reached.")

    return nu


@jit
def danby(k, r0, v0, tof, numiter=20, rtol=1e-8):
    """Kepler solver for both elliptic and parabolic orbits based on Danby's
    algorithm.

    Parameters
    ----------
    k : float
        Standard gravitational parameter of the attractor.
    r0 : numpy.ndarray
        Position vector.
    v0 : numpy.ndarray
        Velocity vector.
    tof : float
        Time of flight.
    numiter : int, optional
        Number of iterations, defaults to 20.
    rtol : float, optional
        Relative error for accuracy of the method, defaults to 1e-8.

    Returns
    -------
    rr : numpy.ndarray
        Final position vector.
    vv : numpy.ndarray
        Final velocity vector.

    Notes
    -----
    This algorithm was developed by Danby in his paper *The solution of Kepler
    Equation* with DOI: https://doi.org/10.1007/BF01686811
    """

    # Solve first for eccentricity and mean anomaly
    p, ecc, inc, raan, argp, nu = rv2coe(k, r0, v0)
    nu = danby_coe(k, p, ecc, inc, raan, argp, nu, tof, numiter, rtol)

    return coe2rv(k, p, ecc, inc, raan, argp, nu)
