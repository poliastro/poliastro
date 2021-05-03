"""Initial orbit determination.

"""
from astropy import units as u

from poliastro.core.iod import vallado as vallado_fast

kms = u.km / u.s


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
    short : bool, optional
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

    v0, v = vallado_fast(k_, r0_, r_, tof_, short, numiter, rtol)

    yield v0 << kms, v << kms
