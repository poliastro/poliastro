"""Ephemerides.

"""

from datetime import datetime, timedelta
from functools import wraps

import numpy as np

from ._octave import octave

__all__ = ['planet_ephem']

EPOCH = datetime(2000, 1, 1, 0, 0)

# TODO: Add the number of the moons too
MERCURY = 1
VENUS = 2
EARTH = 3
MARS = 4
JUPITER = 5
SATURN = 6
URANUS = 7
NEPTUNE = 8
PLUTO = 9
SUN = 10

_cache = {}

def _memoize(key, cache):
    def _decorator(func):
        @wraps(func)
        def _inner(jd, body):
            cache_key = "{:s}:{:.1f}:{:d}".format(key, jd, body)
            if cache_key in cache:
                return cache[cache_key]
            ret = func(jd, body)
            _cache[cache_key] = ret
            return ret
        return _inner
    return _decorator

@_memoize("ephem", _cache)
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

    TODO: Refactor to use constants instead of codes
    E.g.: ephem.ephem(jd, ephem.EARTH)

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


def time_from_epoch(dd=None, epoch=EPOCH):
    """Returns timedelta object from datetime subtracting EPOCH.

    If called without arguments, use datetime.today().
    Use module EPOCH by default, but custom epoch may be provided.

    """
    if dd:
        delta = dd - epoch
    else:
        delta = datetime.today() - epoch
    return delta


def delta2date(td, epoch=EPOCH):
    """Returns date object from timedelta or seconds, adding EPOCH.

    TODO: The naming in this module is quite inconsistent.
    TODO: Refactor

    """
    try:
        dd = EPOCH + td
    except TypeError:
        td = timedelta(seconds=td)
        dd = EPOCH + td
    return dd
