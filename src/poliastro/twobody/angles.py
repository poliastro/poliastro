"""Angles and anomalies.

"""
import numpy as np

from astropy import units as u
from poliastro.jit import jit


@jit
def newton(func, x0, fprime, args=(), tol=1.48e-08, maxiter=50):
    if tol <= 0:
        raise ValueError("tol too small (%g <= 0)" % tol)
    if maxiter < 1:
        raise ValueError("maxiter must be greater than 0")
    # Multiply by 1.0 to convert to floating point.  We don't use float(x0)
    # so it still works if x0 is complex.
    p0 = 1.0 * x0
    for iter in range(maxiter):
        fder = fprime(p0, *args)
        fval = func(p0, *args)

        newton_step = fval / fder
        p = p0 - newton_step
        if abs(p - p0) < tol:
            return p
        p0 = p
    msg = "Failed to converge after %d iterations, value is %s" % (maxiter, p)
    raise RuntimeError(msg)


@jit
def _kepler_equation(E, M, ecc):
    return E - ecc * np.sin(E) - M


@jit
def _kepler_equation_prime(E, M, ecc):
    return 1 - ecc * np.cos(E)


@jit
def _kepler_equation_hyper(F, M, ecc):
    return -F + ecc * np.sinh(F) - M


@jit
def _kepler_equation_prime_hyper(F, M, ecc):
    return ecc * np.cosh(F) - 1


@jit
def _kepler_equation_parabolic(D, M, ecc):
    return M_parabolic(ecc, D) - M


@jit
def _kepler_equation_prime_parabolic(D, M, ecc):
    return M_parabolic_prime(ecc, D)


@jit
def M_parabolic(ecc, D, tolerance=1e-16):
    """Computes the Kepler equation r.h.s. in near-parabolic regime

    Parameters
    ----------
    D : float
        Eccentric anomaly (rad).
    ecc : float
        Eccentricity,
    tolerance : float (optional)
        smallness of the last term in series
    Returns
    -------
    M_parabolic : float
        kepler equation r.h.s.

    Notes
    -----
    Taken from Farnocchia, Davide, Davide Bracali Cioci, and Andrea Milani.
    "Robust resolution of Kepler’s equation in all eccentricity regimes."
    Celestial Mechanics and Dynamical Astronomy 116, no. 1 (2013): 21-34.
    """
    x = (ecc - 1.0) / (ecc + 1.0) * (D ** 2)
    small_term = False
    S = 0.0
    k = 0
    while not small_term:
        term = (ecc - 1.0 / (2.0 * k + 3.0)) * (x ** k)
        small_term = np.abs(term) < tolerance
        S += term
        k += 1
    return np.sqrt(2.0 / (1.0 + ecc)) * D + np.sqrt(2.0 / (1.0 + ecc) ** 3) * (D ** 3) * S


@jit
def M_parabolic_prime(ecc, D, tolerance=1e-16):
    """Computes derivative of the Kepler equation r.h.s. in near-parabolic regime

    Parameters
    ----------
    D : float
        Eccentric anomaly (rad).
    ecc : float
        Eccentricity,
    tolerance : float (optional)
        smallness of the last term in series
    Returns
    -------
    M_parabolic : float
        derivative of kepler equation r.h.s.

    Notes
    -----
    Taken from Farnocchia, Davide, Davide Bracali Cioci, and Andrea Milani.
    "Robust resolution of Kepler’s equation in all eccentricity regimes."
    Celestial Mechanics and Dynamical Astronomy 116, no. 1 (2013): 21-34.
    """
    x = (ecc - 1.0) / (ecc + 1.0) * (D ** 2)
    small_term = False
    S_prime = 0.0
    k = 0
    while not small_term:
        term = (ecc - 1.0 / (2.0 * k + 3.0)) * (2 * k + 3.0) * (x ** k)
        small_term = np.abs(term) < tolerance
        S_prime += term
        k += 1
    return np.sqrt(2.0 / (1.0 + ecc)) + np.sqrt(2.0 / (1.0 + ecc) ** 3) * (D ** 2) * S_prime


def D_to_nu(D, ecc):
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
    return 2.0 * np.arctan(D)


def nu_to_D(nu, ecc):
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
    return np.tan(nu / 2.0)


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
    with u.set_enabled_equivalencies(u.dimensionless_angles()):
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
        E = newton(_kepler_equation, M, _kepler_equation_prime,
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
        F = newton(_kepler_equation_hyper, np.arcsinh(M / ecc), _kepler_equation_prime_hyper,
                            args=(M, ecc), maxiter=100)
    return F


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
    with u.set_enabled_equivalencies(u.dimensionless_angles()):
        B = 3.0 * M / 2.0
        A = (B + (1.0 + B ** 2) ** (0.5)) ** (2.0 / 3.0)
        guess = 2 * A * B / (1 + A + A ** 2)
        D = newton(_kepler_equation_parabolic, guess, _kepler_equation_prime_parabolic,
                            args=(M, ecc), maxiter=100)
    return D


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
    with u.set_enabled_equivalencies(u.dimensionless_angles()):
        M = _kepler_equation_parabolic(D, 0.0 * u.rad, ecc)
    return M


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
    if ecc > 1 + delta:
        F = M_to_F(M, ecc)
        nu = F_to_nu(F, ecc)
    elif ecc < 1 - delta:
        E = M_to_E(M, ecc)
        nu = E_to_nu(E, ecc)
    else:
        D = M_to_D(M, ecc)
        nu = D_to_nu(D, ecc)
    return nu


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
    if ecc > 1 + delta:
        F = nu_to_F(nu, ecc)
        M = F_to_M(F, ecc)
    elif ecc < 1 - delta:
        E = nu_to_E(nu, ecc)
        M = E_to_M(E, ecc)
    else:
        D = nu_to_D(nu, ecc)
        M = D_to_M(D, ecc)
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
