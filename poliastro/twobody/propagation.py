# coding: utf-8
"""Propagation algorithms.

"""
import numpy as np

from poliastro.jit import jit
from poliastro.util import dot
from poliastro.stumpff import c2, c3


def func_twobody(t0, u, a, k, kwargs):
    """Differential equation for generic twobody motion.

    Parameters
    ----------
    u : ndarray
        Six component state vector [x, y, z, vx, vy, vz] (km, km/s).
    t0 : float
        Time.
    a : function(t0, u, k, **kwargs)
        Non Keplerian acceleration (km/s2).
    k : float
        Standard gravitational parameter.
    kwargs : dict
        Dictionary of extra parameters for both the objective function and
        the user-supplied acceleration function.

    """
    ax, ay, az = a(t0, u, k, **kwargs)

    x, y, z, vx, vy, vz = u
    r3 = (x**2 + y**2 + z**2)**1.5

    du = np.array([
        vx,
        vy,
        vz,
        -k * x / r3 + ax,
        -k * y / r3 + ay,
        -k * z / r3 + az
    ])
    return du


def kepler(k, r0, v0, tof, numiter=35, rtol=1e-10):
    """Propagates Keplerian orbit.

    Parameters
    ----------
    k : float
        Gravitational constant of main attractor (km^3 / s^2).
    r0 : array
        Initial position (km).
    v0 : array
        Initial velocity (km).
    tof : float
        Time of flight (s).
    numiter : int, optional
        Maximum number of iterations, default to 35.
    rtol : float, optional
        Maximum relative error permitted, default to 1e-10.

    Raises
    ------
    RuntimeError
        If the algorithm didn't converge.

    Notes
    -----
    This algorithm is based on Vallado implementation, and does basic Newton
    iteration on the Kepler equation written using universal variables. Battin
    claims his algorithm uses the same amount of memory but is between 40 %
    and 85 % faster.

    """
    # Compute Lagrange coefficients
    f, g, fdot, gdot = _kepler(k, r0, v0, tof, numiter, rtol)

    assert np.abs(f * gdot - fdot * g - 1) < 1e-5  # Fixed tolerance

    # Return position and velocity vectors
    r = f * r0 + g * v0
    v = fdot * r0 + gdot * v0

    return r, v


@jit
def _kepler(k, r0, v0, tof, numiter, rtol):
    # Cache some results
    dot_r0v0 = dot(r0, v0)
    norm_r0 = dot(r0, r0) ** .5
    sqrt_mu = k**.5
    alpha = -dot(v0, v0) / k + 2 / norm_r0

    # First guess
    if alpha > 0:
        # Elliptic orbit
        xi_new = sqrt_mu * tof * alpha
    elif alpha < 0:
        # Hyperbolic orbit
        xi_new = (np.sign(tof) * (-1 / alpha)**.5 *
                  np.log((-2 * k * alpha * tof) /
                         (dot_r0v0 + np.sign(tof) *
                          np.sqrt(-k / alpha) * (1 - norm_r0 * alpha))))
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
        norm_r = (xi * xi * c2_psi +
                  dot_r0v0 / sqrt_mu * xi * (1 - psi * c3_psi) +
                  norm_r0 * (1 - psi * c2_psi))
        xi_new = xi + (sqrt_mu * tof - xi * xi * xi * c3_psi -
                       dot_r0v0 / sqrt_mu * xi * xi * c2_psi -
                       norm_r0 * xi * (1 - psi * c3_psi)) / norm_r
        if abs(np.divide(xi_new - xi, xi_new)) < rtol or abs(xi_new - xi) < rtol:
            break
        else:
            count += 1
    else:
        raise RuntimeError("Maximum number of iterations reached")

    # Compute Lagrange coefficients
    f = 1 - xi**2 / norm_r0 * c2_psi
    g = tof - xi**3 / sqrt_mu * c3_psi

    gdot = 1 - xi**2 / norm_r * c2_psi
    fdot = sqrt_mu / (norm_r * norm_r0) * xi * (psi * c3_psi - 1)

    return f, g, fdot, gdot
