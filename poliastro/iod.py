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
    # TODO: Raise custom exceptions
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


def target2(k, r0, v0, r0_tg, v0_tg, tof):
    """Solves the targetting problem.

    """
    v1t, v2t, dv1, dv2, error = _astiod.target(r0, v0, r0_tg, v0_tg, 'S', 'K', tof, k)
    error = error.strip().decode('ascii')
    if error != 'ok':
        raise RuntimeError("There was an error: {}".format(error))
    return v1t, v2t, dv1, dv2, error


def target(k, r0, v0, r0_tg, v0_tg, tof):
    # Propagate target trajectory
    rf, vf, err = twobody.kepler(k, r0_tg, v0_tg, tof)
    # Solve Lambert's problem
    va_s, vb_s = lambert(k, r0, rf, 'short', tof)
    va_l, vb_l = lambert(k, r0, rf, 'long', tof)
    # Compute delta velocity for short and long trajectories
    dva_s = va_s - v0
    dvb_s = vf - vb_s
    dv_s = norm(dva_s) + norm(dvb_s)
    dva_l = va_l - v0
    dvb_l = vf - vb_l
    dv_l = norm(dva_l) + norm(dvb_l)
    # Short or long solution
    if dv_s < dv_l:
        va, vb = va_s, vb_s
        deltav = dv_s
    else:
        va, vb = va_l, vb_l
        deltav = dv_l

    return rf, vf, va, vb, deltav
