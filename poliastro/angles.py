"""Angles and anomalies.

"""

import numpy as np

from ._ast2body import newtonm as M2nu

__all__ = ['fp_angle', 'M2nu', 'nu2M']

M2nu.__doc__ = """Eccentric and true anomaly from mean anomaly.

Parameters
----------
ecc : float
    Eccentricity.
M : float
    Mean anomaly (rad).

Returns
-------
e0 : float
    Eccentric anomaly (rad).
nu : float
    True anomaly (rad).

Examples
--------
>>> _, nu = M2nu(0.06, np.radians(30.0))
>>> np.degrees(nu)
33.673284930211658

"""


def nu2M(ecc, nu):
    """Eccentric and mean anomaly from true anomaly.

    Parameters
    ----------
    ecc : float
        Eccentricity.
    nu : float
        True anomaly (rad).

    Returns
    -------
    e0 : float
        Eccentric anomaly (rad).
    M : float
        Mean anomaly (rad).

    """
    e0 = nu2E(ecc, nu)
    M = e0 - ecc * np.sin(e0)
    return e0, M


def nu2E(ecc, nu):
    e0 = 2 * np.arctan(np.sqrt((1 - ecc) / (1 + ecc)) * np.tan(nu / 2))
    return e0


def fp_angle(ecc, nu):
    """Computes flight path angle.

    Parameters
    ----------
    ecc : float
        Eccentricity.
    nu : float
        True anomaly (rad).

    From Vallado 2007, pp. 113.

    """
    return np.arctan2(ecc * np.sin(nu), 1 + ecc * np.cos(nu))
