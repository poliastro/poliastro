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
    return 2.0 * np.arctan(D)


@jit
def nu_to_D(nu):
    return np.tan(nu / 2.0)


@jit
def nu_to_E(nu, ecc):
    E = 2 * np.arctan(np.sqrt((1 - ecc) / (1 + ecc)) * np.tan(nu / 2))
    return E


@jit
def nu_to_F(nu, ecc):
    F = np.log(
        (np.sqrt(ecc + 1) + np.sqrt(ecc - 1) * np.tan(nu / 2))
        / (np.sqrt(ecc + 1) - np.sqrt(ecc - 1) * np.tan(nu / 2))
    )
    return F


@jit
def E_to_nu(E, ecc):
    nu = 2 * np.arctan(np.sqrt((1 + ecc) / (1 - ecc)) * np.tan(E / 2))
    return nu


@jit
def F_to_nu(F, ecc):
    nu = 2 * np.arctan(
        (np.exp(F) * np.sqrt(ecc + 1) - np.sqrt(ecc + 1))
        / (np.exp(F) * np.sqrt(ecc - 1) + np.sqrt(ecc - 1))
    )
    return nu


@jit
def M_to_E(M, ecc):
    E = newton("elliptic", M, args=(M, ecc))
    return E


@jit
def M_to_F(M, ecc):
    F = newton("hyperbolic", np.arcsinh(M / ecc), args=(M, ecc), maxiter=100)
    return F


@jit
def M_to_D(M, ecc):
    B = 3.0 * M / 2.0
    A = (B + (1.0 + B ** 2) ** 0.5) ** (2.0 / 3.0)
    guess = 2 * A * B / (1 + A + A ** 2)
    D = newton("parabolic", guess, args=(M, ecc), maxiter=100)
    return D


@jit
def E_to_M(E, ecc):
    M = _kepler_equation(E, 0.0, ecc)
    return M


@jit
def F_to_M(F, ecc):
    M = _kepler_equation_hyper(F, 0.0, ecc)
    return M


@jit
def D_to_M(D, ecc):
    M = _kepler_equation_parabolic(D, 0.0, ecc)
    return M


@jit
def M_to_nu(M, ecc, delta=1e-2):
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
    return np.arctan2(ecc * np.sin(nu), 1 + ecc * np.cos(nu))
