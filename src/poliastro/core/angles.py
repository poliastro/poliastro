import numpy as np

from ._jit import jit


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
    x = (ecc - 1.0) / (ecc + 1.0) * (D ** 2)
    small_term = False
    S = 0.0
    k = 0
    while not small_term:
        term = (ecc - 1.0 / (2.0 * k + 3.0)) * (x ** k)
        small_term = np.abs(term) < tolerance
        S += term
        k += 1
    return (
        np.sqrt(2.0 / (1.0 + ecc)) * D + np.sqrt(2.0 / (1.0 + ecc) ** 3) * (D ** 3) * S
    )


@jit
def M_parabolic_prime(ecc, D, tolerance=1e-16):
    x = (ecc - 1.0) / (ecc + 1.0) * (D ** 2)
    small_term = False
    S_prime = 0.0
    k = 0
    while not small_term:
        term = (ecc - 1.0 / (2.0 * k + 3.0)) * (2 * k + 3.0) * (x ** k)
        small_term = np.abs(term) < tolerance
        S_prime += term
        k += 1
    return (
        np.sqrt(2.0 / (1.0 + ecc))
        + np.sqrt(2.0 / (1.0 + ecc) ** 3) * (D ** 2) * S_prime
    )


@jit
def newton(regime, x0, args=(), tol=1.48e-08, maxiter=50):
    p0 = 1.0 * x0
    for iter in range(maxiter):
        if regime == "parabolic":
            fval = _kepler_equation_parabolic(p0, *args)
            fder = _kepler_equation_prime_parabolic(p0, *args)
        elif regime == "hyperbolic":
            fval = _kepler_equation_hyper(p0, *args)
            fder = _kepler_equation_prime_hyper(p0, *args)
        else:
            fval = _kepler_equation(p0, *args)
            fder = _kepler_equation_prime(p0, *args)

        newton_step = fval / fder
        p = p0 - newton_step
        if abs(p - p0) < tol:
            return p
        p0 = p
    return 1.0


@jit
def D_to_nu(D):
    r"""True anomaly from parabolic eccentric anomaly.

    .. math::

        \nu = 2 \cdot \arctan{(D)}

    Parameters
    ----------
    D : ~astropy.units.Quantity
        Eccentric anomaly.

    Returns
    -------
    nu : ~astropy.units.Quantity
        True anomaly.

    Note
    ----
    Taken from Farnocchia, Davide, Davide Bracali Cioci, and Andrea Milani.
    "Robust resolution of Kepler’s equation in all eccentricity regimes."
    Celes
    """

    return 2.0 * np.arctan(D)


@jit
def nu_to_D(nu):
    r"""Parabolic eccentric anomaly from true anomaly.

    .. math::
        D = \tan{\frac{\nu}{2}}

    Parameters
    ----------
    nu : ~astropy.units.Quantity
        True anomaly.

    Returns
    -------
    D : ~astropy.units.Quantity
        Hyperbolic eccentric anomaly.

    Note
    ----
    Taken from Farnocchia, Davide, Davide Bracali Cioci, and Andrea Milani.
    "Robust resolution of Kepler’s equation in all eccentricity regimes."
    Celestial Mechanics and Dynamical Astronomy 116, no. 1 (2013): 21-34.
    """
    return np.tan(nu / 2.0)


@jit
def nu_to_E(nu, ecc):
    r"""Eccentric anomaly from true anomaly.

    .. versionadded:: 0.4.0

    .. math::
        E = 2\arctan{\sqrt{\frac{1-e}{1+e}}\tan{\frac{\nu}{2}}}

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

    beta = ecc / (1 + np.sqrt(1 - (ecc ** 2)))
    E = nu - 2 * np.arctan(beta * np.sin(nu) / (1 + beta * np.cos(nu)))
    return E


@jit
def nu_to_F(nu, ecc):
    r"""Hyperbolic eccentric anomaly from true anomaly.

    .. math::
        F = ln{\left ( \frac{\sin{(\nu)}\sqrt{e^{2}-1} + \cos{\nu} + e}{1+e\cos{(\nu)}} \right )}

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

    Note
    -----
    Taken from Curtis, H. (2013). *Orbital mechanics for engineering students*. 167

    """
    F = np.log(
        (np.sqrt(ecc + 1) + np.sqrt(ecc - 1) * np.tan(nu / 2))
        / (np.sqrt(ecc + 1) - np.sqrt(ecc - 1) * np.tan(nu / 2))
    )
    return F


@jit
def E_to_nu(E, ecc):
    r"""True anomaly from eccentric anomaly.

    .. versionadded:: 0.4.0

    .. math::
        \nu = 2\arctan{\left ( \sqrt{\frac{1+e}{1-e}}\tan{\frac{E}{2}} \right )}

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
    beta = ecc / (1 + np.sqrt((1 - ecc) * (1 + ecc)))
    nu = E + 2 * np.arctan(beta * np.sin(E) / (1 - beta * np.cos(E)))
    return nu


@jit
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
    nu = 2 * np.arctan(
        (np.exp(F) * np.sqrt(ecc + 1) - np.sqrt(ecc + 1))
        / (np.exp(F) * np.sqrt(ecc - 1) + np.sqrt(ecc - 1))
    )
    return nu


@jit
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
    E = newton("elliptic", M, args=(M, ecc))
    return E


@jit
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
    F = newton("hyperbolic", np.arcsinh(M / ecc), args=(M, ecc), maxiter=100)
    return F


@jit
def M_to_D(M, ecc):
    """Parabolic eccentric anomaly from mean anomaly.

    Parameters
    ----------
    M : ~astropy.units.Quantity
        Mean anomaly.
    ecc : ~astropy.units.Quantity
        Eccentricity (>1).

    Returns
    -------
    D : ~astropy.units.Quantity
        Parabolic eccentric anomaly.

    """
    B = 3.0 * M / 2.0
    A = (B + (1.0 + B ** 2) ** 0.5) ** (2.0 / 3.0)
    guess = 2 * A * B / (1 + A + A ** 2)
    D = newton("parabolic", guess, args=(M, ecc), maxiter=100)
    return D


@jit
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
    M = _kepler_equation(E, 0.0, ecc)
    return M


@jit
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
    M = _kepler_equation_hyper(F, 0.0, ecc)
    return M


@jit
def D_to_M(D, ecc):
    """Mean anomaly from eccentric anomaly.

    Parameters
    ----------
    D : ~astropy.units.Quantity
        Parabolic eccentric anomaly.
    ecc : ~astropy.units.Quantity
        Eccentricity.

    Returns
    -------
    M : ~astropy.units.Quantity
        Mean anomaly.

    """
    M = _kepler_equation_parabolic(D, 0.0, ecc)
    return M


@jit
def M_to_nu(M, ecc, delta=1e-2):
    """True anomaly from mean anomaly.

    .. versionadded:: 0.4.0

    Parameters
    ----------
    M : ~astropy.units.Quantity
        Mean anomaly.
    ecc : ~astropy.units.Quantity
        Eccentricity.
    delta : float (optional)
        threshold of near-parabolic regime definition (from Davide Farnocchia et al)
    Returns
    -------
    nu : ~astropy.units.Quantity
        True anomaly.

    Examples
    --------
    >>> import astropy.units as u
    >>> M_to_nu(30.0 * u.deg, 0.06 * u.one)
    -1.5355907568869296

    """
    if ecc > 1 + delta:
        F = M_to_F(M, ecc)
        nu = F_to_nu(F, ecc)
    elif ecc < 1 - delta:
        E = M_to_E(M, ecc)
        nu = E_to_nu(E, ecc)
    else:
        D = M_to_D(M, ecc)
        nu = D_to_nu(D)
    return nu


@jit
def nu_to_M(nu, ecc, delta=1e-2):
    """Mean anomaly from true anomaly.

    .. versionadded:: 0.4.0

    Parameters
    ----------
    nu : ~astropy.units.Quantity
        True anomaly.
    ecc : ~astropy.units.Quantity
        Eccentricity.
    delta : float (optional)
        threshold of near-parabolic regime definition (from Davide Farnocchia et al)

    Returns
    -------
    M : ~astropy.units.Quantity
        Mean anomaly.

    """
    if ecc > 1 + delta:
        F = nu_to_F(nu, ecc)
        M = F_to_M(F, ecc)
    elif ecc < 1 - delta:
        E = nu_to_E(nu, ecc)
        M = E_to_M(E, ecc)
    else:
        D = nu_to_D(nu)
        M = D_to_M(D, ecc)
    return M


@jit
def fp_angle(nu, ecc):
    r"""Returns the flight path angle.

    .. math::
        \gamma = \arctan{\frac{e\sin{\theta}}{1 + e\cos{\theta}}}

    Parameters
    ----------
    nu: ~astropy.units.Quantity
        True anomaly
    ecc: ~astropy.units.Quantity
        Eccentricity

    Returns
    fp_angle: ~astropy.units.Quantity
        Flight path angle

    Note
    -----
    Algorithm taken from Vallado 2007, pp. 113.

    """
    return np.arctan2(ecc * np.sin(nu), 1 + ecc * np.cos(nu))
