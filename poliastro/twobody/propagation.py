# coding: utf-8
"""Propagation algorithms.

"""

import numpy as np

from . import _ast2body


def kepler(k, r0, v0, tof):
    """Propagates orbit.

    This is a wrapper around kepler from ast2body.for.

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

    Raises
    ------
    RuntimeError
        If the status of the subroutine is not 'ok'.

    """
    r0 = np.asarray(r0).astype(np.float)
    v0 = np.asarray(v0).astype(np.float)
    tof = float(tof)
    assert r0.shape == (3,)
    assert v0.shape == (3,)
    r, v, error = _ast2body.kepler(r0, v0, tof, k)
    error = error.strip().decode('ascii')
    if error != 'ok':
        raise RuntimeError("There was an error: {}".format(error))
    return r, v
