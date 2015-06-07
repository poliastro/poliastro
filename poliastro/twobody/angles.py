# coding: utf-8
"""Angles and anomalies.

"""

import numpy as np

from scipy import optimize


def _kepler_equation(E, M, ecc):
    return E - ecc * np.sin(E) - M


def _kepler_equation_prime(E, M, ecc):
    return 1 - ecc * np.cos(E)


def nu_to_E(nu, ecc):
    """Eccentric anomaly from true anomaly.

    .. versionadded:: 0.4.0

    Parameters
    ----------
    nu : float
        True anomaly (rad).
    ecc : float
        Eccentricity.

    Returns
    -------
    E : float
        Eccentric anomaly.

    """
    E = 2 * np.arctan(np.sqrt((1 - ecc) / (1 + ecc)) * np.tan(nu / 2))
    return E


def E_to_nu(E, ecc):
    """True anomaly from eccentric anomaly.

    .. versionadded:: 0.4.0

    Parameters
    ----------
    E : float
        Eccentric anomaly (rad).
    ecc : float
        Eccentricity.

    Returns
    -------
    nu : float
        True anomaly (rad).

    """
    nu = 2 * np.arctan(np.sqrt((1 + ecc) / (1 - ecc)) * np.tan(E / 2))
    return nu


def M_to_E(M, ecc):
    """Eccentric anomaly from true anomaly.

    .. versionadded:: 0.4.0

    Parameters
    ----------
    M : float
        Mean anomaly (rad).
    ecc : float
        Eccentricity.

    Returns
    -------
    E : float
        Eccentric anomaly.

    """
    E = optimize.newton(_kepler_equation, M, _kepler_equation_prime,
                        args=(M, ecc))
    return E


def E_to_M(E, ecc):
    """Mean anomaly from eccentric anomaly.

    .. versionadded:: 0.4.0

    Parameters
    ----------
    E : float
        Eccentric anomaly (rad).
    ecc : float
        Eccentricity.

    Returns
    -------
    M : float
        Mean anomaly (rad).

    """
    M = _kepler_equation(E, 0.0, ecc)
    return M


def M_to_nu(M, ecc):
    """True anomaly from mean anomaly.

    .. versionadded:: 0.4.0

    Parameters
    ----------
    M : float
        Mean anomaly (rad).
    ecc : float
        Eccentricity.

    Returns
    -------
    nu : float
        True anomaly (rad).

    Examples
    --------
    >>> nu = M2nu(np.radians(30.0), 0.06)
    >>> np.rad2deg(nu)
    33.673284930211658

    """
    E = M_to_E(M, ecc)
    nu = E_to_nu(E, ecc)
    return nu


def nu_to_M(nu, ecc):
    """Mean anomaly from true anomaly.

    .. versionadded:: 0.4.0

    Parameters
    ----------
    nu : float
        True anomaly (rad).
    ecc : float
        Eccentricity.

    Returns
    -------
    M : float
        Mean anomaly (rad).

    """
    E = nu_to_E(nu, ecc)
    M = E_to_M(E, ecc)
    return M


def fp_angle(nu, ecc):
    """Flight path angle.

    .. versionadded:: 0.4.0

    Parameters
    ----------
    nu : float
        True anomaly (rad).
    ecc : float
        Eccentricity.

    Notes
    -----
    Algorithm taken from Vallado 2007, pp. 113.

    """
    return np.arctan2(ecc * np.sin(nu), 1 + ecc * np.cos(nu))
