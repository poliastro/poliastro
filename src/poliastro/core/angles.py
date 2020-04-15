import numpy as np

from ._jit import jit


@jit
def _kepler_equation(E, M, ecc):
    return E_to_M(E, ecc) - M


@jit
def _kepler_equation_prime(E, M, ecc):
    return 1 - ecc * np.cos(E)


@jit
def _kepler_equation_hyper(F, M, ecc):
    return F_to_M(F, ecc) - M


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
def S_x(ecc, x, atol=1e-12, maxiter=100):
    S = 0
    k = 0
    while k < maxiter:
        S_old = S
        S += (ecc - 1 / (2 * k + 3)) * x ** k
        k += 1
        if abs(S - S_old) < atol:
            return S
    else:
        raise RuntimeError("Function did not converge")


@jit
def dS_x_alt(ecc, x, atol=1e-12, maxiter=100):
    # Notice that this is not exactly
    # the partial derivative of S with respect to D,
    # but the result of arranging the terms
    # in section 4.2 of Farnocchia et al. 2013
    S = 0
    k = 0
    while k < maxiter:
        S_old = S
        S += (ecc - 1 / (2 * k + 3)) * (2 * k + 3) * x ** k
        k += 1
        if abs(S - S_old) < atol:
            return S
    else:
        raise RuntimeError("Function did not converge")


@jit
def d2S_x_alt(ecc, x, atol=1e-12, maxiter=100):
    # Notice that this is not exactly
    # the second partial derivative of S with respect to D,
    # but the result of arranging the terms
    # in section 4.2 of Farnocchia et al. 2013
    # Also, notice that we are not using this function yet
    S = 0
    k = 0
    while k < maxiter:
        S_old = S
        S += (ecc - 1 / (2 * k + 3)) * (2 * k + 3) * (2 * k + 2) * x ** k
        k += 1
        if abs(S - S_old) < atol:
            return S
    else:
        raise RuntimeError("Function did not converge")


@jit
def M_parabolic(ecc, D, tolerance=1e-16):
    x = (ecc - 1.0) / (ecc + 1.0) * (D ** 2)
    assert abs(x) < 1
    S = S_x(ecc, x, tolerance)
    return (
        np.sqrt(2.0 / (1.0 + ecc)) * D + np.sqrt(2.0 / (1.0 + ecc) ** 3) * (D ** 3) * S
    )


@jit
def M_parabolic_prime(ecc, D, tolerance=1e-16):
    x = (ecc - 1.0) / (ecc + 1.0) * (D ** 2)
    assert abs(x) < 1
    S = dS_x_alt(ecc, x, tolerance)
    return np.sqrt(2.0 / (1.0 + ecc)) + np.sqrt(2.0 / (1.0 + ecc) ** 3) * (D ** 2) * S


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
    D : float
        Eccentric anomaly.

    Returns
    -------
    nu : float
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
    nu : float
        True anomaly in radians.

    Returns
    -------
    D : float
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
    nu : float
        True anomaly in radians.
    ecc : float
        Eccentricity.

    Returns
    -------
    E : float
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
    nu : float
        True anomaly in radians.
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
    E : float
        Eccentric anomaly in radians.
    ecc : float
        Eccentricity.

    Returns
    -------
    nu : float
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
    F : float
        Hyperbolic eccentric anomaly.
    ecc : float
        Eccentricity (>1).

    Returns
    -------
    nu : float
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
    M : float
        Mean anomaly in radians.
    ecc : float
        Eccentricity.

    Returns
    -------
    E : float
        Eccentric anomaly.

    """
    E0 = M
    E = newton("elliptic", E0, args=(M, ecc))
    return E


@jit
def M_to_F(M, ecc):
    """Hyperbolic eccentric anomaly from mean anomaly.

    Parameters
    ----------
    M : float
        Mean anomaly in radians.
    ecc : float
        Eccentricity (>1).

    Returns
    -------
    F : float
        Hyperbolic eccentric anomaly.

    """
    F0 = np.arcsinh(M / ecc)
    F = newton("hyperbolic", F0, args=(M, ecc), maxiter=100)
    return F


@jit
def M_to_D(M):
    """Parabolic eccentric anomaly from mean anomaly.

    Parameters
    ----------
    M : float
        Mean anomaly in radians.

    Returns
    -------
    D : float
        Parabolic eccentric anomaly.

    Notes
    -----
    This uses the analytical solution of Barker's equation,
    see Battin, 1987.

    """
    B = 3.0 * M / 2.0
    A = (B + (1.0 + B ** 2) ** 0.5) ** (2.0 / 3.0)
    D = 2 * A * B / (1 + A + A ** 2)
    return D


@jit
def M_to_D_near_parabolic(M, ecc):
    """Parabolic eccentric anomaly from mean anomaly, near parabolic case.

    Parameters
    ----------
    M : float
        Mean anomaly in radians.
    ecc : float
        Eccentricity (~1).

    Returns
    -------
    D : float
        Parabolic eccentric anomaly.

    """
    D0 = M_to_D(M)
    D = newton("parabolic", D0, args=(M, ecc), maxiter=100)
    return D


@jit
def E_to_M(E, ecc):
    """Mean anomaly from eccentric anomaly.

    .. versionadded:: 0.4.0

    Parameters
    ----------
    E : float
        Eccentric anomaly in radians.
    ecc : float
        Eccentricity.

    Returns
    -------
    M : float
        Mean anomaly.

    """
    M = E - ecc * np.sin(E)
    return M


@jit
def F_to_M(F, ecc):
    """Mean anomaly from eccentric anomaly.

    Parameters
    ----------
    F : float
        Hyperbolic eccentric anomaly.
    ecc : float
        Eccentricity (>1).

    Returns
    -------
    M : float
        Mean anomaly.

    """
    M = ecc * np.sinh(F) - F
    return M


@jit
def D_to_M(D, ecc):
    """Mean anomaly from eccentric anomaly.

    Parameters
    ----------
    D : float
        Parabolic eccentric anomaly.
    ecc : float
        Eccentricity.

    Returns
    -------
    M : float
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
    M : float
        Mean anomaly in radians.
    ecc : float
        Eccentricity.
    delta : float (optional)
        threshold of near-parabolic regime definition (from Davide Farnocchia et al)

    Returns
    -------
    nu : float
        True anomaly.

    Examples
    --------
    >>> from numpy import radians, degrees
    >>> degrees(M_to_nu(radians(30.0), 0.06))
    33.67328493021166

    """
    if ecc > 1 + delta:
        F = M_to_F(M, ecc)
        nu = F_to_nu(F, ecc)
    elif ecc < 1 - delta:
        E = M_to_E(M, ecc)
        nu = E_to_nu(E, ecc)
    else:
        D = M_to_D_near_parabolic(M, ecc)
        nu = D_to_nu(D)
    return nu


@jit
def nu_to_M(nu, ecc, delta=1e-2):
    """Mean anomaly from true anomaly.

    .. versionadded:: 0.4.0

    Parameters
    ----------
    nu : float
        True anomaly in radians.
    ecc : float
        Eccentricity.
    delta : float (optional)
        Threshold of near-parabolic regime definition (from Davide Farnocchia et al)

    Returns
    -------
    M : float
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
    nu: float
        True anomaly in radians.
    ecc: float
        Eccentricity.

    Returns
    fp_angle: float
        Flight path angle

    Note
    -----
    Algorithm taken from Vallado 2007, pp. 113.

    """
    return np.arctan2(ecc * np.sin(nu), 1 + ecc * np.cos(nu))
