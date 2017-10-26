"""Angles and anomalies.

"""
import numpy as np

from astropy import units as u

from scipy import optimize


def _kepler_equation(E, M, ecc):
    return E - ecc * np.sin(E) - M


def _kepler_equation_prime(E, M, ecc):
    return 1 - ecc * np.cos(E)


def _kepler_equation_hyper(F, M, ecc):
    return -F + ecc * np.sinh(F) - M


def _kepler_equation_prime_hyper(F, M, ecc):
    return ecc * np.cosh(F) - 1


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


def nu_to_F(nu, ecc):
    """Hyperbolic eccentric anomaly from true anomaly.

    Parameters
    ----------
    nu : float
        True anomaly (rad).
    ecc : float
        Eccentricity (>1).

    Returns
    -------
    F : float
        Hyperbolic eccentric anomaly.

    Note
    -----
    Taken from Curtis, H. (2013). *Orbital mechanics for engineering students*. 167

    """
    F = np.log((np.sqrt(ecc + 1) + np.sqrt(ecc - 1) * np.tan(nu / 2)) /
               (np.sqrt(ecc + 1) - np.sqrt(ecc - 1) * np.tan(nu / 2))) * u.rad
    return F


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


def F_to_nu(F, ecc):
    """True anomaly from hyperbolic eccentric anomaly.

    Parameters
    ----------
    F : float
        Hyperbolic eccentric anomaly (rad).
    ecc : float
        Eccentricity (>1).

    Returns
    -------
    nu : float
        True anomaly (rad).

    """
    nu = 2 * np.arctan((np.exp(F) * np.sqrt(ecc + 1) - np.sqrt(ecc + 1)) /
                       (np.exp(F) * np.sqrt(ecc - 1) + np.sqrt(ecc - 1)))
    return nu


def M_to_E(M, ecc):
    """Eccentric anomaly from mean anomaly.

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
    with u.set_enabled_equivalencies(u.dimensionless_angles()):
        E = optimize.newton(_kepler_equation, M, _kepler_equation_prime,
                            args=(M, ecc))
    return E


def M_to_F(M, ecc):
    """Hyperbolic eccentric anomaly from mean anomaly.

    Parameters
    ----------
    M : float
        Mean anomaly (rad).
    ecc : float
        Eccentricity (>1).

    Returns
    -------
    F : float
        Hyperbolic eccentric anomaly.

    """
    with u.set_enabled_equivalencies(u.dimensionless_angles()):
        F = optimize.newton(_kepler_equation_hyper, M, _kepler_equation_prime_hyper,
                            args=(M, ecc))
    return F


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
    with u.set_enabled_equivalencies(u.dimensionless_angles()):
        M = _kepler_equation(E, 0.0 * u.rad, ecc)
    return M


def F_to_M(F, ecc):
    """Mean anomaly from eccentric anomaly.

    Parameters
    ----------
    F : float
        Hyperbolic eccentric anomaly (rad).
    ecc : float
        Eccentricity (>1).

    Returns
    -------
    M : float
        Mean anomaly (rad).

    """
    with u.set_enabled_equivalencies(u.dimensionless_angles()):
        M = _kepler_equation_hyper(F, 0.0 * u.rad, ecc)
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
    >>> nu = M_to_nu(np.radians(30.0), 0.06)
    >>> np.rad2deg(nu)
    33.673284930211658

    """
    if ecc > 1:
        F = M_to_F(M, ecc)
        nu = F_to_nu(F, ecc)
    else:
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
    if ecc > 1:
        F = nu_to_F(nu, ecc)
        M = F_to_M(F, ecc)
    else:
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

    Note
    -----
    Algorithm taken from Vallado 2007, pp. 113.

    """
    return np.arctan2(ecc * np.sin(nu), 1 + ecc * np.cos(nu))
