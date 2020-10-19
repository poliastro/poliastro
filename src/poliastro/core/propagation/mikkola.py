import numpy as np

from .._jit import jit
from ..angles import D_to_nu, E_to_M, E_to_nu, F_to_M, F_to_nu, nu_to_E, nu_to_F
from ..elements import coe2rv, rv2coe


@jit
def mikkola_coe(k, p, ecc, inc, raan, argp, nu, tof):

    a = p / (1 - ecc ** 2)
    n = np.sqrt(k / np.abs(a) ** 3)

    # Solve for specific geometrical case
    if ecc < 1.0:
        # Equation (9a)
        alpha = (1 - ecc) / (4 * ecc + 1 / 2)
        M0 = E_to_M(nu_to_E(nu, ecc), ecc)
    else:
        alpha = (ecc - 1) / (4 * ecc + 1 / 2)
        M0 = F_to_M(nu_to_F(nu, ecc), ecc)

    M = M0 + n * tof
    beta = M / 2 / (4 * ecc + 1 / 2)

    # Equation (9b)
    if beta >= 0:
        z = (beta + np.sqrt(beta ** 2 + alpha ** 3)) ** (1 / 3)
    else:
        z = (beta - np.sqrt(beta ** 2 + alpha ** 3)) ** (1 / 3)

    s = z - alpha / z

    # Apply initial correction
    if ecc < 1.0:
        ds = -0.078 * s ** 5 / (1 + ecc)
    else:
        ds = 0.071 * s ** 5 / (1 + 0.45 * s ** 2) / (1 + 4 * s ** 2) / ecc

    s += ds

    # Solving for the true anomaly
    if ecc < 1.0:
        E = M + ecc * (3 * s - 4 * s ** 3)
        f = E - ecc * np.sin(E) - M
        f1 = 1.0 - ecc * np.cos(E)
        f2 = ecc * np.sin(E)
        f3 = ecc * np.cos(E)
        f4 = -f2
        f5 = -f3
    else:
        E = 3 * np.log(s + np.sqrt(1 + s ** 2))
        f = -E + ecc * np.sinh(E) - M
        f1 = -1.0 + ecc * np.cosh(E)
        f2 = ecc * np.sinh(E)
        f3 = ecc * np.cosh(E)
        f4 = f2
        f5 = f3

    # Apply Taylor expansion
    u1 = -f / f1
    u2 = -f / (f1 + 0.5 * f2 * u1)
    u3 = -f / (f1 + 0.5 * f2 * u2 + (1.0 / 6.0) * f3 * u2 ** 2)
    u4 = -f / (
        f1 + 0.5 * f2 * u3 + (1.0 / 6.0) * f3 * u3 ** 2 + (1.0 / 24.0) * f4 * (u3 ** 3)
    )
    u5 = -f / (
        f1
        + f2 * u4 / 2
        + f3 * (u4 * u4) / 6.0
        + f4 * (u4 * u4 * u4) / 24.0
        + f5 * (u4 * u4 * u4 * u4) / 120.0
    )

    E += u5

    if ecc < 1.0:
        nu = E_to_nu(E, ecc)
    else:
        if ecc == 1.0:
            # Parabolic
            nu = D_to_nu(E)
        else:
            # Hyperbolic
            nu = F_to_nu(E, ecc)

    return nu


@jit
def mikkola(k, r0, v0, tof, rtol=None):
    """Raw algorithm for Mikkola's Kepler solver.

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
    nu = mikkola_coe(k, p, ecc, inc, raan, argp, nu, tof)

    return coe2rv(k, p, ecc, inc, raan, argp, nu)


@jit
def mikkola_tofs_coe(k, p, ecc, inc, raan, argp, nu, tofs):
    nus = np.zeros_like(tofs)
    for i, tof in tofs.ndenumerate():
        nus[i] = mikkola_coe(k, p, ecc, inc, raan, argp, nu, tof)
    return nus


@jit
def mikkola_tofs(k, r0, v0, tofs, rtol=None):
    p, ecc, inc, raan, argp, nu = rv2coe(k, r0, v0)
    return mikkola_tofs_coe(k, p, ecc, inc, raan, argp, nu, tofs)
