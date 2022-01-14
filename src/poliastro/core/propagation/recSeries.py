import numpy as np
from numba import njit as jit

from poliastro.core.angles import E_to_M, nu_to_E, E_to_nu
from poliastro.core.elements import coe2rv, rv2coe


@jit
def recSeries_coe(k, p, ecc, inc, raan, argp, nu, tof, numiter=20, order=8):

    semi_axis_a = p / (1 - ecc ** 2)
    n = np.sqrt(k / np.abs(semi_axis_a) ** 3)

    
    if ecc == 0:
        # Solving for circular orbit
        
        M0 = nu
        M = M0 + n * tof
        nu = M - 2 * np.pi * np.floor(M / 2 / np.pi)
        return nu

    elif ecc < 1.0:
        # Solving for elliptical orbit
        
        M0 = E_to_M(nu_to_E(nu, ecc), ecc)
        M = M0 + n * tof
        M = M - 2 * np.pi * np.floor(M / 2 / np.pi)
        
        E = M
        for i in range(0,order):
            E = M + ecc*np.sin(E)
        
        return E_to_nu(E,ecc)
        
        
    else:
        raise ValueError("Parabolic/Hyperbolic orbits not supported.")

    return nu


@jit
def recSeries(k, r0, v0, tof, numiter=20, order=8):
    """Kepler solver for elliptical orbits with recursive series approximation
    method. The order of the series is a user defined parameter.

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
    order : int, optional
        Order of recursion, defaults to 8.

    Returns
    -------
    rr : numpy.ndarray
        Final position vector.
    vv : numpy.ndarray
        Final velocity vector.

    Notes
    -----
    This algorithm uses series developed in the paper *Recursive solution to
    Kepler’s problem for elliptical orbits - application in robust 
    Newton-Raphson and co-planar closest approach estimation* 
    with DOI: http://dx.doi.org/10.13140/RG.2.2.18578.58563/1
    """

    # Solve first for eccentricity and mean anomaly
    p, ecc, inc, raan, argp, nu = rv2coe(k, r0, v0)
    nu = recSeries_coe(k, p, ecc, inc, raan, argp, nu, tof, numiter, order)

    return coe2rv(k, p, ecc, inc, raan, argp, nu)
