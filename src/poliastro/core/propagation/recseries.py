import numpy as np
from numba import njit as jit

from poliastro.core.angles import E_to_M, E_to_nu, nu_to_E
from poliastro.core.elements import coe2rv, rv2coe


@jit
def recseries_coe(
    k,
    p,
    ecc,
    inc,
    raan,
    argp,
    nu,
    tof,
    method="rtol",
    order=8,
    numiter=100,
    rtol=1e-8,
):

    # semi-major axis
    semi_axis_a = p / (1 - ecc**2)
    # mean angular motion
    n = np.sqrt(k / np.abs(semi_axis_a) ** 3)

    if ecc == 0:
        # Solving for circular orbit

        # compute initial mean anoamly
        M0 = nu  # For circular orbit (M = E = nu)
        # final mean anaomaly
        M = M0 + n * tof
        # snapping anomaly to [0,pi] range
        nu = M - 2 * np.pi * np.floor(M / 2 / np.pi)

        return nu

    elif ecc < 1.0:
        # Solving for elliptical orbit

        # compute initial mean anoamly
        M0 = E_to_M(nu_to_E(nu, ecc), ecc)
        # final mean anaomaly
        M = M0 + n * tof
        # snapping anomaly to [0,pi] range
        M = M - 2 * np.pi * np.floor(M / 2 / np.pi)

        # set recursion iteration
        if method == "rtol":
            Niter = numiter
        elif method == "order":
            Niter = order
        else:
            raise ValueError(
                "Unknown recursion termination method ('rtol','order')."
            )

        # compute eccentric anomaly through recursive series
        E = M + ecc  # Using initial guess from vallado to improve convergence
        for i in range(0, Niter):
            En = M + ecc * np.sin(E)
            # check for break condition
            if method == "rtol" and (abs(En - E) / abs(E)) < rtol:
                break
            E = En

        return E_to_nu(E, ecc)

    else:
        # Parabolic/Hyperbolic orbits are not supported
        raise ValueError("Parabolic/Hyperbolic orbits not supported.")

    return nu


@jit
def recseries(k, r0, v0, tof, method="rtol", order=8, numiter=100, rtol=1e-8):
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
    method : str
        Type of termination method ('rtol','order')
    order : int, optional
        Order of recursion, defaults to 8.
    numiter : int, optional
        Number of iterations, defaults to 100.
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
    This algorithm uses series discussed in the paper *Recursive solution to
    Kepler’s problem for elliptical orbits - application in robust
    Newton-Raphson and co-planar closest approach estimation*
    with DOI: http://dx.doi.org/10.13140/RG.2.2.18578.58563/1
    """

    # Solve first for eccentricity and mean anomaly
    p, ecc, inc, raan, argp, nu = rv2coe(k, r0, v0)
    nu = recseries_coe(
        k, p, ecc, inc, raan, argp, nu, tof, method, order, numiter, rtol
    )

    return coe2rv(k, p, ecc, inc, raan, argp, nu)
