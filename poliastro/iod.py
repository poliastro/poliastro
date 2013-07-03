"""Initial orbit determination.

"""

import numpy as np
from numpy.linalg import norm

from poliastro import twobody
from . import _astiod

__all__ = ['lambert', 'target']


def lambert(k, r0, rf, tof, short=True):
    """
    Solves the Lambert problem.

    This is a wrapper around LAMBERTUNIV from ASTIOD.FOR.

    Parameters
    ----------
    k : float
        Gravitational constant of main attractor (km^3/s^2).
    r0 : array
        Initial position (km).
    r : array
        Final position (km).
    tof : float
        Time of flight (s).
    short : boolean
        Find out the short path. Default to True.

    Raises
    ------
    RuntimeError
        If the status of the subroutine is not 'ok'.

    Note
    ----
    Multiple revolutions not supported.

    """
    r0 = np.asanyarray(r0).astype(np.float)
    rf = np.asanyarray(rf).astype(np.float)
    tof = float(tof)
    assert r0.shape == (3,)
    assert rf.shape == (3,)
    assert type(short) is bool
    if short:
        dm = "SHORT"
    else:
        dm = "LONG"
    va, vb, error = _astiod.lambertuniv(r0, rf, dm, 'N', tof, k)
    error = error.strip().decode('ascii')
    if error != 'ok':
        raise RuntimeError("There was an error: {}".format(error))
    return va, vb


def target(k, r0, v0, r0_tg, v0_tg, tof):
    # TODO: Docstring
    rf, vf = twobody.kepler(k, r0_tg, v0_tg, tof)
    va_short, vb_short = lambert(k, r0, rf, tof, True)
    va_long, vb_long = lambert(k, r0, rf, tof, False)
    dva_short = va_short - v0
    dvb_short = vf - vb_short
    dv_short = norm(dva_short) + norm(dvb_short)
    dva_long = va_long - v0
    dvb_long = vf - vb_long
    dv_long = norm(dva_long) + norm(dvb_long)
    # Choose short or long solution
    if dv_short < dv_long:
        va, vb = va_short, vb_short
        dva, dvb = dva_short, dvb_short
        #deltav = dv_short
    else:
        va, vb = va_long, vb_long
        dva, dvb = dva_long, dvb_long
        #deltav = dv_long

    return rf, vf, va, vb, dva, dvb
