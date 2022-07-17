import numpy as np
from numba import njit as jit

from poliastro.core.angles import (
    D_to_M,
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
from poliastro.core.elements import coe2rv, rv2coe


@jit
def _kepler_equation_near_parabolic(D, M, ecc):
    return D_to_M_near_parabolic(D, ecc) - M


@jit
def _kepler_equation_prime_near_parabolic(D, M, ecc):
    x = (ecc - 1.0) / (ecc + 1.0) * (D**2)
    assert abs(x) < 1
    S = dS_x_alt(ecc, x)
    return (
        np.sqrt(2.0 / (1.0 + ecc))
        + np.sqrt(2.0 / (1.0 + ecc) ** 3) * (D**2) * S
    )


@jit
def S_x(ecc, x, atol=1e-12):
    assert abs(x) < 1
    S = 0
    k = 0
    while True:
        S_old = S
        S += (ecc - 1 / (2 * k + 3)) * x**k
        k += 1
        if abs(S - S_old) < atol:
            return S


@jit
def dS_x_alt(ecc, x, atol=1e-12):
    # Notice that this is not exactly
    # the partial derivative of S with respect to D,
    # but the result of arranging the terms
    # in section 4.2 of Farnocchia et al. 2013
    assert abs(x) < 1
    S = 0
    k = 0
    while True:
        S_old = S
        S += (ecc - 1 / (2 * k + 3)) * (2 * k + 3) * x**k
        k += 1
        if abs(S - S_old) < atol:
            return S


@jit
def d2S_x_alt(ecc, x, atol=1e-12):
    # Notice that this is not exactly
    # the second partial derivative of S with respect to D,
    # but the result of arranging the terms
    # in section 4.2 of Farnocchia et al. 2013
    # Also, notice that we are not using this function yet
    assert abs(x) < 1
    S = 0
    k = 0
    while True:
        S_old = S
        S += (ecc - 1 / (2 * k + 3)) * (2 * k + 3) * (2 * k + 2) * x**k
        k += 1
        if abs(S - S_old) < atol:
            return S


@jit
def D_to_M_near_parabolic(D, ecc):
    x = (ecc - 1.0) / (ecc + 1.0) * (D**2)
    assert abs(x) < 1
    S = S_x(ecc, x)
    return (
        np.sqrt(2.0 / (1.0 + ecc)) * D
        + np.sqrt(2.0 / (1.0 + ecc) ** 3) * (D**3) * S
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
    tol : float, optional
        Absolute tolerance for Newton convergence.
    maxiter : int, optional
        Maximum number of iterations for Newton convergence.

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
def delta_t_from_nu(nu, ecc, k=1.0, q=1.0, delta=1e-2):
    """Time elapsed since periapsis for given true anomaly.

    Parameters
    ----------
    nu : float
        True anomaly.
    ecc : float
        Eccentricity.
    k : float
        Gravitational parameter.
    q : float
        Periapsis distance.
    delta : float
        Parameter that controls the size of the near parabolic region.

    Returns
    -------
    delta_t : float
        Time elapsed since periapsis.

    """
    assert -np.pi <= nu < np.pi
    if ecc < 1 - delta:
        # Strong elliptic
        E = nu_to_E(nu, ecc)  # (-pi, pi]
        M = E_to_M(E, ecc)  # (-pi, pi]
        n = np.sqrt(k * (1 - ecc) ** 3 / q**3)
    elif 1 - delta <= ecc < 1:
        E = nu_to_E(nu, ecc)  # (-pi, pi]
        if delta <= 1 - ecc * np.cos(E):
            # Strong elliptic
            M = E_to_M(E, ecc)  # (-pi, pi]
            n = np.sqrt(k * (1 - ecc) ** 3 / q**3)
        else:
            # Near parabolic
            D = nu_to_D(nu)  # (-∞, ∞)
            # If |nu| is far from pi this result is bounded
            # because the near parabolic region shrinks in its vicinity,
            # otherwise the eccentricity is very close to 1
            # and we are really far away
            M = D_to_M_near_parabolic(D, ecc)
            n = np.sqrt(k / (2 * q**3))
    elif ecc == 1:
        # Parabolic
        D = nu_to_D(nu)  # (-∞, ∞)
        M = D_to_M(D)  # (-∞, ∞)
        n = np.sqrt(k / (2 * q**3))
    elif 1 + ecc * np.cos(nu) < 0:
        # Unfeasible region
        return np.nan
    elif 1 < ecc <= 1 + delta:
        # NOTE: Do we need to wrap nu here?
        # For hyperbolic orbits, it should anyway be in
        # (-arccos(-1 / ecc), +arccos(-1 / ecc))
        F = nu_to_F(nu, ecc)  # (-∞, ∞)
        if delta <= ecc * np.cosh(F) - 1:
            # Strong hyperbolic
            M = F_to_M(F, ecc)  # (-∞, ∞)
            n = np.sqrt(k * (ecc - 1) ** 3 / q**3)
        else:
            # Near parabolic
            D = nu_to_D(nu)  # (-∞, ∞)
            M = D_to_M_near_parabolic(D, ecc)  # (-∞, ∞)
            n = np.sqrt(k / (2 * q**3))
    elif 1 + delta < ecc:
        # Strong hyperbolic
        F = nu_to_F(nu, ecc)  # (-∞, ∞)
        M = F_to_M(F, ecc)  # (-∞, ∞)
        n = np.sqrt(k * (ecc - 1) ** 3 / q**3)
    else:
        raise RuntimeError

    return M / n


@jit
def nu_from_delta_t(delta_t, ecc, k=1.0, q=1.0, delta=1e-2):
    """True anomaly for given elapsed time since periapsis.

    Parameters
    ----------
    delta_t : float
        Time elapsed since periapsis.
    ecc : float
        Eccentricity.
    k : float
        Gravitational parameter.
    q : float
        Periapsis distance.
    delta : float
        Parameter that controls the size of the near parabolic region.

    Returns
    -------
    nu : float
        True anomaly.

    """
    if ecc < 1 - delta:
        # Strong elliptic
        n = np.sqrt(k * (1 - ecc) ** 3 / q**3)
        M = n * delta_t
        # This might represent several revolutions,
        # so we wrap the true anomaly
        E = M_to_E((M + np.pi) % (2 * np.pi) - np.pi, ecc)
        nu = E_to_nu(E, ecc)
    elif 1 - delta <= ecc < 1:
        E_delta = np.arccos((1 - delta) / ecc)
        # We compute M assuming we are in the strong elliptic case
        # and verify later
        n = np.sqrt(k * (1 - ecc) ** 3 / q**3)
        M = n * delta_t
        # We check against abs(M) because E_delta could also be negative
        if E_to_M(E_delta, ecc) <= abs(M):
            # Strong elliptic, proceed
            # This might represent several revolutions,
            # so we wrap the true anomaly
            E = M_to_E((M + np.pi) % (2 * np.pi) - np.pi, ecc)
            nu = E_to_nu(E, ecc)
        else:
            # Near parabolic, recompute M
            n = np.sqrt(k / (2 * q**3))
            M = n * delta_t
            D = M_to_D_near_parabolic(M, ecc)
            nu = D_to_nu(D)
    elif ecc == 1:
        # Parabolic
        n = np.sqrt(k / (2 * q**3))
        M = n * delta_t
        D = M_to_D(M)
        nu = D_to_nu(D)
    elif 1 < ecc <= 1 + delta:
        F_delta = np.arccosh((1 + delta) / ecc)
        # We compute M assuming we are in the strong hyperbolic case
        # and verify later
        n = np.sqrt(k * (ecc - 1) ** 3 / q**3)
        M = n * delta_t
        # We check against abs(M) because F_delta could also be negative
        if F_to_M(F_delta, ecc) <= abs(M):
            # Strong hyperbolic, proceed
            F = M_to_F(M, ecc)
            nu = F_to_nu(F, ecc)
        else:
            # Near parabolic, recompute M
            n = np.sqrt(k / (2 * q**3))
            M = n * delta_t
            D = M_to_D_near_parabolic(M, ecc)
            nu = D_to_nu(D)
    # elif 1 + delta < ecc:
    else:
        # Strong hyperbolic
        n = np.sqrt(k * (ecc - 1) ** 3 / q**3)
        M = n * delta_t
        F = M_to_F(M, ecc)
        nu = F_to_nu(F, ecc)

    return nu


@jit
def farnocchia_coe(k, p, ecc, inc, raan, argp, nu, tof):

    q = p / (1 + ecc)

    delta_t0 = delta_t_from_nu(nu, ecc, k, q)
    delta_t = delta_t0 + tof

    return nu_from_delta_t(delta_t, ecc, k, q)


@jit
def farnocchia_rv(k, r0, v0, tof):
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
    r0 : numpy.ndarray
        Initial position vector wrt attractor center.
    v0 : numpy.ndarray
        Initial velocity vector.
    tof : float
        Time of flight (s).

    Notes
    -----
    This method takes initial :math:`\vec{r}, \vec{v}`, calculates classical orbit parameters,
    increases mean anomaly and performs inverse transformation to get final :math:`\vec{r}, \vec{v}`
    The logic is based on formulae (4), (6) and (7) from http://dx.doi.org/10.1007/s10569-013-9476-9

    """

    # get the initial true anomaly and orbit parameters that are constant over time
    p, ecc, inc, raan, argp, nu0 = rv2coe(k, r0, v0)
    nu = farnocchia_coe(k, p, ecc, inc, raan, argp, nu0, tof)

    return coe2rv(k, p, ecc, inc, raan, argp, nu)
