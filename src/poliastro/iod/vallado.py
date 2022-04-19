"""Initial orbit determination.

"""
from astropy import units as u

from poliastro.core.iod import vallado as vallado_fast

kms = u.km / u.s


def lambert(
    k, r0, r, tof, M=0, prograde=True, lowpath=True, numiter=35, rtol=1e-8
):
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
    M : int, optional
        Number of full revolutions, default to 0.
    prograde: boolean
        Controls the desired inclination of the transfer orbit.
    lowpath: boolean
        If `True` or `False`, gets the transfer orbit whose vacant focus is
        below or above the chord line, respectively.
    numiter : int, optional
        Maximum number of iterations, default to 35.
    rtol : float, optional
        Relative tolerance of the algorithm, default to 1e-8.

    Returns
    -------
    v0, v : tuple
        Pair of velocity solutions.

    Raises
    ------
    RuntimeError
        If it was not possible to compute the orbit.

    Notes
    -----
    This uses the universal variable approach found in Battin, Mueller & White
    with the bisection iteration suggested by Vallado. Multiple revolutions
    not supported.

    """
    k_ = k.to_value(u.km**3 / u.s**2)
    r0_ = r0.to_value(u.km)
    r_ = r.to_value(u.km)
    tof_ = tof.to_value(u.s)

    v0, v = vallado_fast(
        k_, r0_, r_, tof_, M, prograde, lowpath, numiter, rtol
    )

    return v0 << kms, v << kms
