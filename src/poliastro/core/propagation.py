import numpy as np

from poliastro.core.angles import M_to_nu, nu_to_M
from poliastro.core.elements import coe2rv, rv2coe
from poliastro.core.stumpff import c2, c3

from ._jit import jit


@jit
def mean_motion(k, r0, v0, tof):
    r"""Propagates orbit using mean motion

    .. versionadded:: 0.9.0

    Parameters
    ----------
    orbit : ~poliastro.twobody.orbit.Orbit
        the Orbit object to propagate.
    tof : float
        Time of flight (s).

    Notes
    -----
    This method takes initial :math:`\vec{r}, \vec{v}`, calculates classical orbit parameters,
    increases mean anomaly and performs inverse transformation to get final :math:`\vec{r}, \vec{v}`
    The logic is based on formulae (4), (6) and (7) from http://dx.doi.org/10.1007/s10569-013-9476-9

    """

    # get the initial true anomaly and orbit parameters that are constant over time
    p, ecc, inc, raan, argp, nu0 = rv2coe(k, r0, v0)

    # get the initial mean anomaly
    M0 = nu_to_M(nu0, ecc)
    # strong elliptic or strong hyperbolic orbits
    if np.abs(ecc - 1.0) > 1e-2:
        a = p / (1.0 - ecc ** 2)
        # given the initial mean anomaly, calculate mean anomaly
        # at the end, mean motion (n) equals sqrt(mu / |a^3|)
        M = M0 + tof * np.sqrt(k / np.abs(a ** 3))
        nu = M_to_nu(M, ecc)

    # near-parabolic orbit
    else:
        q = p * np.abs(1.0 - ecc) / np.abs(1.0 - ecc ** 2)
        # mean motion n = sqrt(mu / 2 q^3) for parabolic orbit
        M = M0 + tof * np.sqrt(k / 2.0 / (q ** 3))
        nu = M_to_nu(M, ecc)

    return coe2rv(k, p, ecc, inc, raan, argp, nu)


@jit
def kepler(k, r0, v0, tof, numiter):
    # Cache some results
    dot_r0v0 = np.dot(r0, v0)
    norm_r0 = np.dot(r0, r0) ** 0.5
    sqrt_mu = k ** 0.5
    alpha = -np.dot(v0, v0) / k + 2 / norm_r0

    # First guess
    if alpha > 0:
        # Elliptic orbit
        xi_new = sqrt_mu * tof * alpha
    elif alpha < 0:
        # Hyperbolic orbit
        xi_new = (
            np.sign(tof)
            * (-1 / alpha) ** 0.5
            * np.log(
                (-2 * k * alpha * tof)
                / (
                    dot_r0v0
                    + np.sign(tof) * np.sqrt(-k / alpha) * (1 - norm_r0 * alpha)
                )
            )
        )
    else:
        # Parabolic orbit
        # (Conservative initial guess)
        xi_new = sqrt_mu * tof / norm_r0

    # Newton-Raphson iteration on the Kepler equation
    count = 0
    while count < numiter:
        xi = xi_new
        psi = xi * xi * alpha
        c2_psi = c2(psi)
        c3_psi = c3(psi)
        norm_r = (
            xi * xi * c2_psi
            + dot_r0v0 / sqrt_mu * xi * (1 - psi * c3_psi)
            + norm_r0 * (1 - psi * c2_psi)
        )
        xi_new = (
            xi
            + (
                sqrt_mu * tof
                - xi * xi * xi * c3_psi
                - dot_r0v0 / sqrt_mu * xi * xi * c2_psi
                - norm_r0 * xi * (1 - psi * c3_psi)
            )
            / norm_r
        )
        if abs(xi_new - xi) < 1e-7:
            break
        else:
            count += 1
    else:
        raise RuntimeError("Maximum number of iterations reached")

    # Compute Lagrange coefficients
    f = 1 - xi ** 2 / norm_r0 * c2_psi
    g = tof - xi ** 3 / sqrt_mu * c3_psi

    gdot = 1 - xi ** 2 / norm_r * c2_psi
    fdot = sqrt_mu / (norm_r * norm_r0) * xi * (psi * c3_psi - 1)

    return f, g, fdot, gdot
