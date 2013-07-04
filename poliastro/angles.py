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

"""


def fp_angle(ecc, nu):
    """Computes flight path angle.

    From Vallado 2007, pp. 113.

    """
    return np.arctan2(ecc * np.sin(nu), 1 + ecc * np.cos(nu))


def nu2M(ecc, nu):
    E0 = nu2E(ecc, nu)
    M = E0 - ecc * np.sin(E0)
    return E0, M
