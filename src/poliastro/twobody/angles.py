"""Angles and anomalies.

"""
from astropy import units as u

from poliastro.core.angles import (
    D_to_M as D_to_M_fast,
    D_to_nu as D_to_nu_fast,
    E_to_M as E_to_M_fast,
    E_to_nu as E_to_nu_fast,
    F_to_M as F_to_M_fast,
    F_to_nu as F_to_nu_fast,
    M_to_D as M_to_D_fast,
    M_to_E as M_to_E_fast,
    M_to_F as M_to_F_fast,
    fp_angle as fp_angle_fast,
    nu_to_D as nu_to_D_fast,
    nu_to_E as nu_to_E_fast,
    nu_to_F as nu_to_F_fast,
)


@u.quantity_input(D=u.rad)
def D_to_nu(D):
    """True anomaly from parabolic eccentric anomaly.

    Parameters
    ----------
    D : ~astropy.units.Quantity
        Eccentric anomaly.

    Returns
    -------
    nu : ~astropy.units.Quantity
        True anomaly.

    Notes
    -----
    Taken from Farnocchia, Davide, Davide Bracali Cioci, and Andrea Milani.
    "Robust resolution of Kepler’s equation in all eccentricity regimes."
    Celestial Mechanics and Dynamical Astronomy 116, no. 1 (2013): 21-34.
    """
    return (D_to_nu_fast(D.to_value(u.rad)) * u.rad).to(D.unit)


@u.quantity_input(nu=u.rad)
def nu_to_D(nu):
    """Parabolic eccentric anomaly from true anomaly.

    Parameters
    ----------
    nu : ~astropy.units.Quantity
        True anomaly.

    Returns
    -------
    D : ~astropy.units.Quantity
        Hyperbolic eccentric anomaly.

    Notes
    -----
    Taken from Farnocchia, Davide, Davide Bracali Cioci, and Andrea Milani.
    "Robust resolution of Kepler’s equation in all eccentricity regimes."
    Celestial Mechanics and Dynamical Astronomy 116, no. 1 (2013): 21-34.
    """
    return (nu_to_D_fast(nu.to_value(u.rad)) * u.rad).to(nu.unit)


@u.quantity_input(nu=u.rad, ecc=u.one)
def nu_to_E(nu, ecc):
    """Eccentric anomaly from true anomaly.

    .. versionadded:: 0.4.0

    Parameters
    ----------
    nu : ~astropy.units.Quantity
        True anomaly.
    ecc : ~astropy.units.Quantity
        Eccentricity.

    Returns
    -------
    E : ~astropy.units.Quantity
        Eccentric anomaly.

    """
    return (nu_to_E_fast(nu.to_value(u.rad), ecc.value) * u.rad).to(nu.unit)


@u.quantity_input(nu=u.rad, ecc=u.one)
def nu_to_F(nu, ecc):
    """Hyperbolic eccentric anomaly from true anomaly.

    Parameters
    ----------
    nu : ~astropy.units.Quantity
        True anomaly.
    ecc : ~astropy.units.Quantity
        Eccentricity (>1).

    Returns
    -------
    F : ~astropy.units.Quantity
        Hyperbolic eccentric anomaly.

    Notes
    -----
    Taken from Curtis, H. (2013). *Orbital mechanics for engineering students*. 167

    """
    return (nu_to_F_fast(nu.to_value(u.rad), ecc.value) * u.rad).to(nu.unit)


@u.quantity_input(E=u.rad, ecc=u.one)
def E_to_nu(E, ecc):
    """True anomaly from eccentric anomaly.

    .. versionadded:: 0.4.0

    Parameters
    ----------
    E : ~astropy.units.Quantity
        Eccentric anomaly.
    ecc : ~astropy.units.Quantity
        Eccentricity.

    Returns
    -------
    nu : ~astropy.units.Quantity
        True anomaly.

    """
    return (E_to_nu_fast(E.to_value(u.rad), ecc.value) * u.rad).to(E.unit)


@u.quantity_input(F=u.rad, ecc=u.one)
def F_to_nu(F, ecc):
    """True anomaly from hyperbolic eccentric anomaly.

    Parameters
    ----------
    F : ~astropy.units.Quantity
        Hyperbolic eccentric anomaly.
    ecc : ~astropy.units.Quantity
        Eccentricity (>1).

    Returns
    -------
    nu : ~astropy.units.Quantity
        True anomaly.

    """
    return (F_to_nu_fast(F.to_value(u.rad), ecc.value) * u.rad).to(F.unit)


@u.quantity_input(M=u.rad, ecc=u.one)
def M_to_E(M, ecc):
    """Eccentric anomaly from mean anomaly.

    .. versionadded:: 0.4.0

    Parameters
    ----------
    M : ~astropy.units.Quantity
        Mean anomaly.
    ecc : ~astropy.units.Quantity
        Eccentricity.

    Returns
    -------
    E : ~astropy.units.Quantity
        Eccentric anomaly.

    """
    return (M_to_E_fast(M.to_value(u.rad), ecc.value) * u.rad).to(M.unit)


@u.quantity_input(M=u.rad, ecc=u.one)
def M_to_F(M, ecc):
    """Hyperbolic eccentric anomaly from mean anomaly.

    Parameters
    ----------
    M : ~astropy.units.Quantity
        Mean anomaly.
    ecc : ~astropy.units.Quantity
        Eccentricity (>1).

    Returns
    -------
    F : ~astropy.units.Quantity
        Hyperbolic eccentric anomaly.

    """
    return (M_to_F_fast(M.to_value(u.rad), ecc.value) * u.rad).to(M.unit)


@u.quantity_input(M=u.rad, ecc=u.one)
def M_to_D(M):
    """Parabolic eccentric anomaly from mean anomaly.

    Parameters
    ----------
    M : ~astropy.units.Quantity
        Mean anomaly.

    Returns
    -------
    D : ~astropy.units.Quantity
        Parabolic eccentric anomaly.

    """
    return (M_to_D_fast(M.to_value(u.rad)) * u.rad).to(M.unit)


@u.quantity_input(E=u.rad, ecc=u.one)
def E_to_M(E, ecc):
    """Mean anomaly from eccentric anomaly.

    .. versionadded:: 0.4.0

    Parameters
    ----------
    E : ~astropy.units.Quantity
        Eccentric anomaly.
    ecc : ~astropy.units.Quantity
        Eccentricity.

    Returns
    -------
    M : ~astropy.units.Quantity
        Mean anomaly.

    """
    return (E_to_M_fast(E.to_value(u.rad), ecc.value) * u.rad).to(E.unit)


@u.quantity_input(F=u.rad, ecc=u.one)
def F_to_M(F, ecc):
    """Mean anomaly from eccentric anomaly.

    Parameters
    ----------
    F : ~astropy.units.Quantity
        Hyperbolic eccentric anomaly.
    ecc : ~astropy.units.Quantity
        Eccentricity (>1).

    Returns
    -------
    M : ~astropy.units.Quantity
        Mean anomaly.

    """
    return (F_to_M_fast(F.to_value(u.rad), ecc.value) * u.rad).to(F.unit)


@u.quantity_input(D=u.rad, ecc=u.one)
def D_to_M(D):
    """Mean anomaly from eccentric anomaly.

    Parameters
    ----------
    D : ~astropy.units.Quantity
        Parabolic eccentric anomaly.

    Returns
    -------
    M : ~astropy.units.Quantity
        Mean anomaly.

    """
    return (D_to_M_fast(D.to_value(u.rad)) * u.rad).to(D.unit)


@u.quantity_input(nu=u.rad, ecc=u.one)
def fp_angle(nu, ecc):
    """Flight path angle.

    .. versionadded:: 0.4.0

    Parameters
    ----------
    nu : ~astropy.units.Quantity
        True anomaly.
    ecc : ~astropy.units.Quantity
        Eccentricity.

    Notes
    -----
    Algorithm taken from Vallado 2007, pp. 113.

    """
    return (fp_angle_fast(nu.to_value(u.rad), ecc.value) * u.rad).to(nu.unit)
