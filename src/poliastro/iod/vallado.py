"""Initial orbit determination.

"""
import numpy as np
from astropy import units as u

from numba import njit as jit

from poliastro.stumpff import c2, c3


def lambert(k, r0, r, tof, short=True, numiter=35, rtol=1e-8):
    """Solves the Lambert problem.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    k : ~astropy.units.Quantity
        Gravitational constant of main attractor (km^3 / s^2).
    r0 : ~astropy.units.Quantity
        Initial position (km).
    r : ~astropy.units.Quantity
        Final position (km).
    tof : ~astropy.units.Quantity
        Time of flight (s).
    short : boolean, optional
        Find out the short path, default to True. If False, find long path.
    numiter : int, optional
        Maximum number of iterations, default to 35.
    rtol : float, optional
        Relative tolerance of the algorithm, default to 1e-8.

    Raises
    ------
    RuntimeError
        If it was not possible to compute the orbit.

    Note
    ----
    This uses the universal variable approach found in Battin, Mueller & White
    with the bisection iteration suggested by Vallado. Multiple revolutions
    not supported.

    """
    k_ = k.to(u.km ** 3 / u.s ** 2).value
    r0_ = r0.to(u.km).value
    r_ = r.to(u.km).value
    tof_ = tof.to(u.s).value

    v0, v = _lambert(k_, r0_, r_, tof_, short, numiter, rtol)

    yield v0 * u.km / u.s, v * u.km / u.s


@jit
def _lambert(k, r0, r, tof, short, numiter, rtol):
    if short:
        t_m = +1
    else:
        t_m = -1

    norm_r0 = np.dot(r0, r0)**.5
    norm_r = np.dot(r, r)**.5
    cos_dnu = np.dot(r0, r) / (norm_r0 * norm_r)

    A = t_m * (norm_r * norm_r0 * (1 + cos_dnu))**.5

    if A == 0.0:
        raise RuntimeError("Cannot compute orbit, phase angle is 180 degrees")

    psi = 0.0
    psi_low = -4 * np.pi
    psi_up = 4 * np.pi

    count = 0
    while count < numiter:
        y = norm_r0 + norm_r + A * (psi * c3(psi) - 1) / c2(psi)**.5
        if A > 0.0 and y < 0.0:
            # Readjust xi_low until y > 0.0
            # Translated directly from Vallado
            while y < 0.0:
                psi_low = psi
                psi = (0.8 * (1.0 / c3(psi)) *
                       (1.0 - (norm_r0 + norm_r) * np.sqrt(c2(psi)) / A))
                y = norm_r0 + norm_r + A * (psi * c3(psi) - 1) / c2(psi)**.5

        xi = np.sqrt(y / c2(psi))
        tof_new = (xi**3 * c3(psi) + A * np.sqrt(y)) / np.sqrt(k)

        # Convergence check
        if np.abs((tof_new - tof) / tof) < rtol:
            break
        else:
            count += 1
            # Bisection check
            if tof_new <= tof:
                psi_low = psi
            else:
                psi_up = psi
            psi = (psi_up + psi_low) / 2
    else:
        raise RuntimeError("Maximum number of iterations reached")

    f = 1 - y / norm_r0
    g = A * np.sqrt(y / k)

    gdot = 1 - y / norm_r

    v0 = (r - f * r0) / g
    v = (gdot * r - r0) / g

    return v0, v
