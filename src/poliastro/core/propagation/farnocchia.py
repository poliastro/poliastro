import numpy as np

from .._jit import jit
from ..angles import (
    D_to_nu,
    E_to_M,
    E_to_nu,
    F_to_M,
    F_to_nu,
    M_to_D,
    M_to_E,
    M_to_F,
    nu_to_D,
    nu_to_E,
    nu_to_F,
)
from ..elements import coe2rv, rv2coe


@jit
def _kepler_equation_near_parabolic(D, M, ecc):
    return D_to_M_near_parabolic(ecc, D) - M


@jit
def _kepler_equation_prime_near_parabolic(D, M, ecc):
    x = (ecc - 1.0) / (ecc + 1.0) * (D ** 2)
    assert abs(x) < 1
    S = dS_x_alt(ecc, x)
    return np.sqrt(2.0 / (1.0 + ecc)) + np.sqrt(2.0 / (1.0 + ecc) ** 3) * (D ** 2) * S


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
def D_to_M_near_parabolic(ecc, D):
    x = (ecc - 1.0) / (ecc + 1.0) * (D ** 2)
    assert abs(x) < 1
    S = S_x(ecc, x)
    return (
        np.sqrt(2.0 / (1.0 + ecc)) * D + np.sqrt(2.0 / (1.0 + ecc) ** 3) * (D ** 3) * S
    )


@jit
def M_to_D_near_parabolic(M, ecc, tol=1.48e-08, maxiter=50):
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

    for _ in range(maxiter):
        fval = _kepler_equation_near_parabolic(D0, M, ecc)
        fder = _kepler_equation_prime_near_parabolic(D0, M, ecc)

        newton_step = fval / fder
        D = D0 - newton_step
        if abs(D - D0) < tol:
            return D

        D0 = D

    return np.nan


@jit
def D_to_M_near_parabolic_alt(D, ecc):
    """Mean anomaly from eccentric anomaly in the near parabolic region.

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
    M = _kepler_equation_near_parabolic(D, 0.0, ecc)
    return M


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
        # FIXME: This should be
        # M = D_to_M_near_parabolic(D, ecc)
        # but the discrimination of near-parabolic regions here
        # is not correct and the function might not converge
        M = D_to_M_near_parabolic_alt(D, ecc)
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
def farnocchia(k, r0, v0, tof):
    r"""Propagates orbit using mean motion.

    This algorithm depends on the geometric shape of the orbit.
    For the case of the strong elliptic or strong hyperbolic orbits:

    ..  math::

        M = M_{0} + \frac{\mu^{2}}{h^{3}}\left ( 1 -e^{2}\right )^{\frac{3}{2}}t

    .. versionadded:: 0.9.0

    Parameters
    ----------
    k : float
        Standar Gravitational parameter
    r0 : ~astropy.units.Quantity
        Initial position vector wrt attractor center.
    v0 : ~astropy.units.Quantity
        Initial velocity vector.
    tof : float
        Time of flight (s).

    Note
    ----
    This method takes initial :math:`\vec{r}, \vec{v}`, calculates classical orbit parameters,
    increases mean anomaly and performs inverse transformation to get final :math:`\vec{r}, \vec{v}`
    The logic is based on formulae (4), (6) and (7) from http://dx.doi.org/10.1007/s10569-013-9476-9

    """

    # get the initial true anomaly and orbit parameters that are constant over time
    p, ecc, inc, raan, argp, nu0 = rv2coe(k, r0, v0)

    # get the initial mean anomaly
    M0 = nu_to_M(nu0, ecc)
    if np.abs(ecc - 1.0) > 1e-2:
        # strong elliptic or strong hyperbolic orbits
        a = p / (1.0 - ecc ** 2)
        n = np.sqrt(k / np.abs(a ** 3))
    else:
        # near-parabolic orbit
        q = p / np.abs(1.0 + ecc)
        n = np.sqrt(k / 2.0 / (q ** 3))

    M = M0 + tof * n
    nu = M_to_nu(M, ecc)

    return coe2rv(k, p, ecc, inc, raan, argp, nu)
