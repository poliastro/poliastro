"""Ephemerides.

"""

import numpy as np

from ._octave import octave

__all__ = ['planet_ephem']


def planet_ephem(jd, body):
    """Returns the ephemerides of a given body given a Julian date.

    This is a wrapper around `uplanet_2013.m`.

    Parameters
    ----------
    jd : float
        Days elapsed from 1/1/2000 h 00.00, precision up to half a day
    body : int
        Body code

    Returns
    -------
    p, ecc, inc, omega, argp, nu

    """
    assert type(body) is int
    jd = float(jd)
    data_up = octave.call("uplanet_2013", jd, body)
    a, ecc, inc, omega, argp, nu = data_up[0]
    p = a * (1 - ecc ** 2)
    inc = np.deg2rad(inc)
    omega = np.deg2rad(omega)
    argp = np.deg2rad(argp)
    nu = np.deg2rad(nu)
    return p, ecc, inc, omega, argp, nu
