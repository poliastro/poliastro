"""Angles and anomalies.

"""
from astropy import units as u
from poliastro.core.angles import (D_to_nu as D_to_nu_fast, nu_to_D as nu_to_D_fast,
                                   E_to_nu as E_to_nu_fast, nu_to_E as nu_to_E_fast,
                                   F_to_nu as F_to_nu_fast, nu_to_F as nu_to_F_fast,
                                   M_to_D as M_to_D_fast, D_to_M as D_to_M_fast,
                                   M_to_E as M_to_E_fast, E_to_M as E_to_M_fast,
                                   M_to_F as M_to_F_fast, F_to_M as F_to_M_fast,
                                   fp_angle as fp_angle_fast,
                                   M_to_nu as M_to_nu_fast, nu_to_M as nu_to_M_fast)


def D_to_nu(D):
    """True anomaly from parabolic eccentric anomaly.

    Parameters
    ----------
    D : float
        Eccentric anomaly (rad).
    ecc : float
        Eccentricity.

    Returns
    -------
    nu : float
        True anomaly (rad).

    Notes
    -----
    Taken from Farnocchia, Davide, Davide Bracali Cioci, and Andrea Milani.
    "Robust resolution of Kepler’s equation in all eccentricity regimes."
    Celestial Mechanics and Dynamical Astronomy 116, no. 1 (2013): 21-34.
    """
    return D_to_nu_fast(D.to(u.rad).value) * u.rad


def nu_to_D(nu):
    """Parabolic eccentric anomaly from true anomaly.

    Parameters
    ----------
    nu : float
        True anomaly (rad).
    ecc : float
        Eccentricity (>1).

    Returns
    -------
    D : float
        Hyperbolic eccentric anomaly.

    Notes
    -----
    Taken from Farnocchia, Davide, Davide Bracali Cioci, and Andrea Milani.
    "Robust resolution of Kepler’s equation in all eccentricity regimes."
    Celestial Mechanics and Dynamical Astronomy 116, no. 1 (2013): 21-34.
    """
    return nu_to_D_fast(nu.to(u.rad).value) * u.rad


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
    return nu_to_E_fast(nu.to(u.rad).value, ecc.value) * u.rad


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
    return nu_to_F_fast(nu.to(u.rad).value, ecc.value) * u.rad


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
    return E_to_nu_fast(E.to(u.rad).value, ecc.value) * u.rad


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
    return F_to_nu_fast(F.to(u.rad).value, ecc.value) * u.rad


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
    return M_to_E_fast(M.to(u.rad).value, ecc.value) * u.rad


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
    return M_to_F_fast(M.to(u.rad).value, ecc.value) * u.rad


def M_to_D(M, ecc):
    """Parabolic eccentric anomaly from mean anomaly.

    Parameters
    ----------
    M : float
        Mean anomaly (rad).
    ecc : float
        Eccentricity (>1).

    Returns
    -------
    D : float
        Parabolic eccentric anomaly.

    """
    return M_to_D_fast(M.to(u.rad).value, ecc.value) * u.rad


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
    return E_to_M_fast(E.to(u.rad).value, ecc.value) * u.rad


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
    return F_to_M_fast(F.to(u.rad).value, ecc.value) * u.rad


def D_to_M(D, ecc):
    """Mean anomaly from eccentric anomaly.

    Parameters
    ----------
    D : float
        Parabolic eccentric anomaly (rad).
    ecc : float
        Eccentricity.

    Returns
    -------
    M : float
        Mean anomaly (rad).

    """
    return D_to_M_fast(D.to(u.rad).value, ecc.value) * u.rad


def M_to_nu(M, ecc, delta=1e-2):
    """True anomaly from mean anomaly.

    .. versionadded:: 0.4.0

    Parameters
    ----------
    M : float
        Mean anomaly (rad).
    ecc : float
        Eccentricity.
    delta : float (optional)
        threshold of near-parabolic regime definition (from Davide Farnocchia et al)
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
    return M_to_nu_fast(M.to(u.rad).value, ecc.value, delta) * u.rad


def nu_to_M(nu, ecc, delta=1e-2):
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
    return nu_to_M_fast(nu.to(u.rad).value, ecc.value, delta) * u.rad


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
    return fp_angle_fast(nu.to(u.rad).value, ecc.value) * u.rad
