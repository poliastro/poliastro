"""Initial orbit determination.

"""

import numpy as np
import numba

from poliastro.util import dot
from poliastro.stumpff import c2, c3


def lambert(k, r0, r, tof, short=True, numiter=35, rtol=1e-8):
    """Solves the Lambert problem.

    Parameters
    ----------
    k : float
        Gravitational constant of main attractor (km^3 / s^2).
    r0 : array
        Initial position (km).
    r : array
        Final position (km).
    tof : float
        Time of flight (s).
    short : boolean
        Find out the short path, default to True. If False, find long path.

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
    f, g, fdot, gdot = _lambert(k, r0, r, tof, short, numiter, rtol)

    v0 = (r - f * r0) / g
    v = (gdot * r - r0) / g

    return v0, v


@numba.njit
def _lambert(k, r0, r, tof, short, numiter, rtol):
    if short:
        t_m = +1
    else:
        t_m = -1

    norm_r0 = dot(r0, r0)**.5
    norm_r = dot(r, r)**.5
    cos_dnu = dot(r0, r) / (norm_r0 * norm_r)

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
                psi = 0.8 * (1.0 / c3(psi)) * (1.0 - (norm_r0 + norm_r) * np.sqrt(c2(psi)) / A)
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

    return f, g, (f * gdot - 1) / g, gdot
