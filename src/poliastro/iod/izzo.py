"""Izzo's algorithm for Lambert's problem

"""
from astropy import units as u

from poliastro.core.iod import izzo as izzo_fast

kms = u.km / u.s


def lambert(k, r0, r, tof, M=0, numiter=35, rtol=1e-8):
    """Solves the Lambert problem using the Izzo algorithm.

    .. versionadded:: 0.5.0

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
    M : int, optional
        Number of full revolutions, default to 0.
    numiter : int, optional
        Maximum number of iterations, default to 35.
    rtol : float, optional
        Relative tolerance of the algorithm, default to 1e-8.

    Yields
    ------
    v0, v : tuple
        Pair of velocity solutions.

    """
    k_ = k.to_value(u.km**3 / u.s**2)
    r0_ = r0.to_value(u.km)
    r_ = r.to_value(u.km)
    tof_ = tof.to_value(u.s)

    sols = izzo_fast(k_, r0_, r_, tof_, M, numiter, rtol)

    for v0, v in sols:
        yield v0 << kms, v << kms
