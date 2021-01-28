import numpy as np
from numba import njit as jit


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
def newton(regime, x0, args=(), tol=1.48e-08, maxiter=50):
    p0 = 1.0 * x0
    for iter in range(maxiter):
        if regime == "hyperbolic":
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

    return np.nan


@jit
def D_to_nu(D):
    r"""True anomaly from parabolic anomaly.

    Parameters
    ----------
    D : float
        Eccentric anomaly.

    Returns
    -------
    nu : float
        True anomaly.

    Notes
    -----
    From [1]_:

    .. math::

        \nu = 2 \arctan{D}

    """

    return 2.0 * np.arctan(D)


@jit
def nu_to_D(nu):
    r"""Parabolic anomaly from true anomaly.

    Parameters
    ----------
    nu : float
        True anomaly in radians.

    Returns
    -------
    D : float
        Parabolic anomaly.

    Warnings
    --------
    The parabolic anomaly will be continuous in (-∞, ∞)
    only if the true anomaly is in (-π, π].
    No validation or wrapping is performed.

    Notes
    -----
    The treatment of the parabolic case is heterogeneous in the literature,
    and that includes the use of an equivalent quantity to the eccentric anomaly:
    [1]_ calls it "parabolic eccentric anomaly" D,
    [2]_ also uses the letter D but calls it just "parabolic anomaly",
    [3]_ uses the letter B citing indirectly [4]_
    (which however calls it "parabolic time argument"),
    and [5]_ does not bother to define it.

    We use this definition:

    .. math::

        B = \tan{\frac{\nu}{2}}

    References
    ----------
    .. [1] Farnocchia, Davide, Davide Bracali Cioci, and Andrea Milani.
       "Robust resolution of Kepler’s equation in all eccentricity regimes."
    .. [2] Bate, Muller, White.
    .. [3] Vallado, David. "Fundamentals of Astrodynamics and Applications",
       2013.
    .. [4] IAU VIth General Assembly, 1938.
    .. [5] Battin, Richard H. "An introduction to the Mathematics and Methods
       of Astrodynamics, Revised Edition", 1999.

    """
    # TODO: Rename to B
    return np.tan(nu / 2.0)


@jit
def nu_to_E(nu, ecc):
    r"""Eccentric anomaly from true anomaly.

    .. versionadded:: 0.4.0

    Parameters
    ----------
    nu : float
        True anomaly in radians.
    ecc : float
        Eccentricity.

    Returns
    -------
    E : float
        Eccentric anomaly, between -π and π radians.

    Warnings
    --------
    The eccentric anomaly will be between -π and π radians,
    no matter the value of the true anomaly.

    Notes
    -----
    The implementation uses the half-angle formula from [3]_:

    .. math::
        E = 2 \arctan \left ( \sqrt{\frac{1 - e}{1 + e}} \tan{\frac{\nu}{2}} \right)
        \in (-\pi, \pi]

    """
    E = 2 * np.arctan(np.sqrt((1 - ecc) / (1 + ecc)) * np.tan(nu / 2))
    return E


@jit
def nu_to_F(nu, ecc):
    r"""Hyperbolic anomaly from true anomaly.

    Parameters
    ----------
    nu : float
        True anomaly in radians.
    ecc : float
        Eccentricity (>1).

    Returns
    -------
    F : float
        Hyperbolic anomaly.

    Warnings
    --------
    The hyperbolic anomaly will be continuous in (-∞, ∞)
    only if the true anomaly is in (-π, π],
    which should happen anyway
    because the true anomaly is limited for hyperbolic orbits.
    No validation or wrapping is performed.

    Notes
    -----
    The implementation uses the half-angle formula from [3]_:

    .. math::
        F = 2 \operatorname{arctanh} \left( \sqrt{\frac{e-1}{e+1}} \tan{\frac{\nu}{2}} \right)

    """
    F = 2 * np.arctanh(np.sqrt((ecc - 1) / (ecc + 1)) * np.tan(nu / 2))
    return F


@jit
def E_to_nu(E, ecc):
    r"""True anomaly from eccentric anomaly.

    .. versionadded:: 0.4.0

    Parameters
    ----------
    E : float
        Eccentric anomaly in radians.
    ecc : float
        Eccentricity.

    Returns
    -------
    nu : float
        True anomaly, between -π and π radians.

    Warnings
    --------
    The true anomaly will be between -π and π radians,
    no matter the value of the eccentric anomaly.

    Notes
    -----
    The implementation uses the half-angle formula from [3]_:

    .. math::
        \nu = 2 \arctan \left( \sqrt{\frac{1 + e}{1 - e}} \tan{\frac{E}{2}} \right)
        \in (-\pi, \pi]

    """
    nu = 2 * np.arctan(np.sqrt((1 + ecc) / (1 - ecc)) * np.tan(E / 2))
    return nu


@jit
def F_to_nu(F, ecc):
    r"""True anomaly from hyperbolic anomaly.

    Parameters
    ----------
    F : float
        Hyperbolic anomaly.
    ecc : float
        Eccentricity (>1).

    Returns
    -------
    nu : float
        True anomaly.

    Notes
    -----
    The implementation uses the half-angle formula from [3]_:

    .. math::
        \nu = 2 \arctan \left( \sqrt{\frac{e + 1}{e - 1}} \tanh{\frac{F}{2}} \right)
        \in (-\pi, \pi]

    """
    nu = 2 * np.arctan(np.sqrt((ecc + 1) / (ecc - 1)) * np.tanh(F / 2))
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

    Notes
    -----
    This uses a Newton iteration on the Kepler equation.

    """
    assert -np.pi <= M <= np.pi
    if ecc < 0.8:
        E0 = M
    else:
        E0 = np.pi * np.sign(M)
    E = newton("elliptic", E0, args=(M, ecc))
    return E


@jit
def M_to_F(M, ecc):
    """Hyperbolic anomaly from mean anomaly.

    Parameters
    ----------
    M : float
        Mean anomaly in radians.
    ecc : float
        Eccentricity (>1).

    Returns
    -------
    F : float
        Hyperbolic anomaly.

    Notes
    -----
    This uses a Newton iteration on the hyperbolic Kepler equation.

    """
    F0 = np.arcsinh(M / ecc)
    F = newton("hyperbolic", F0, args=(M, ecc), maxiter=100)
    return F


@jit
def M_to_D(M):
    """Parabolic anomaly from mean anomaly.

    Parameters
    ----------
    M : float
        Mean anomaly in radians.

    Returns
    -------
    D : float
        Parabolic anomaly.

    Notes
    -----
    This uses the analytical solution of Barker's equation from [5]_.

    """
    B = 3.0 * M / 2.0
    A = (B + (1.0 + B ** 2) ** 0.5) ** (2.0 / 3.0)
    D = 2 * A * B / (1 + A + A ** 2)
    return D


@jit
def E_to_M(E, ecc):
    r"""Mean anomaly from eccentric anomaly.

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

    Warnings
    --------
    The mean anomaly will be outside of (-π, π]
    if the eccentric anomaly is.
    No validation or wrapping is performed.

    Notes
    -----
    The implementation uses the plain original Kepler equation:

    .. math::
        M = E - e \sin{E}

    """
    M = E - ecc * np.sin(E)
    return M


@jit
def F_to_M(F, ecc):
    r"""Mean anomaly from eccentric anomaly.

    Parameters
    ----------
    F : float
        Hyperbolic anomaly.
    ecc : float
        Eccentricity (>1).

    Returns
    -------
    M : float
        Mean anomaly.

    Notes
    -----
    As noted in [5]_, by manipulating
    the parametric equations of the hyperbola
    we can derive a quantity that is equivalent
    to the eccentric anomaly in the elliptic case:

    .. math::

        M = e \sinh{F} - F

    """
    M = ecc * np.sinh(F) - F
    return M


@jit
def D_to_M(D):
    r"""Mean anomaly from parabolic anomaly.

    Parameters
    ----------
    D : float
        Parabolic anomaly.

    Returns
    -------
    M : float
        Mean anomaly.

    Notes
    -----
    We use this definition:

    .. math::

        M = B + \frac{B^3}{3}

    Notice that M < ν until ν ~ 100 degrees,
    then it reaches π when ν ~ 120 degrees,
    and grows without bounds after that.
    Therefore, it can hardly be called an "anomaly"
    since it is by no means an angle.

    """
    M = D + D ** 3 / 3
    return M


@jit
def fp_angle(nu, ecc):
    r"""Returns the flight path angle.

    Parameters
    ----------
    nu: float
        True anomaly in radians.
    ecc: float
        Eccentricity.

    Returns
    -------
    fp_angle: float
        Flight path angle

    Notes
    -----
    From [3]_, pp. 113:

    .. math::

        \phi = \arctan(\frac {e \sin{\nu}}{1 + e \cos{\nu}})

    """
    return np.arctan2(ecc * np.sin(nu), 1 + ecc * np.cos(nu))
