import numpy as np
from numba import njit as jit

from poliastro.core.angles import E_to_M, E_to_nu, nu_to_E
from poliastro.core.elements import coe2rv, rv2coe


@jit
def gooding_coe(k, p, ecc, inc, raan, argp, nu, tof, numiter=150, rtol=1e-8):
    # TODO: parabolic and hyperbolic not implemented cases
    if ecc >= 1.0:
        raise NotImplementedError(
            "Parabolic/Hyperbolic cases still not implemented in gooding."
        )

    M0 = E_to_M(nu_to_E(nu, ecc), ecc)
    semi_axis_a = p / (1 - ecc**2)
    n = np.sqrt(k / np.abs(semi_axis_a) ** 3)
    M = M0 + n * tof

    # Start the computation
    n = 0
    c = ecc * np.cos(M)
    s = ecc * np.sin(M)
    psi = s / np.sqrt(1 - 2 * c + ecc**2)
    f = 1.0
    while f**2 >= rtol and n <= numiter:
        xi = np.cos(psi)
        eta = np.sin(psi)
        fd = (1 - c * xi) + s * eta
        fdd = c * eta + s * xi
        f = psi - fdd
        psi = psi - f * fd / (fd**2 - 0.5 * f * fdd)
        n += 1

    E = M + psi
    return E_to_nu(E, ecc)


@jit
def gooding(k, r0, v0, tof, numiter=150, rtol=1e-8):
    """Solves the Elliptic Kepler Equation with a cubic convergence and
    accuracy better than 10e-12 rad is normally achieved. It is not valid for
    eccentricities equal or higher than 1.0.

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
        Number of iterations, defaults to 150.
    rtol : float, optional
        Relative error for accuracy of the method, defaults to 1e-8.

    Returns
    -------
    rr : numpy.ndarray
        Final position vector.
    vv : numpy.ndarray
        Final velocity vector.
    Note
    ----
    Original paper for the algorithm: https://doi.org/10.1007/BF01238923
    """

    # Solve first for eccentricity and mean anomaly
    p, ecc, inc, raan, argp, nu = rv2coe(k, r0, v0)
    nu = gooding_coe(k, p, ecc, inc, raan, argp, nu, tof, numiter, rtol)

    return coe2rv(k, p, ecc, inc, raan, argp, nu)
