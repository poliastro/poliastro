"""Izzo's algorithm for Lambert's problem

"""
import numpy as np
from numpy import pi
from astropy import units as u

from scipy.optimize import newton  # Halley iteration

# from poliastro.jit import jit
from poliastro.util import norm

try:
    from hyper.hyper import hyp2f1
except ImportError:
    from scipy.special import hyp2f1


def lambert(k, r0, r, tof, M=0, numiter=35, rtol=1e-8):
    """Solves the Lambert problem using the Izzo algorithm.

    .. versionadded:: 0.5.0

    Parameters
    ----------
    k : Quantity
        Gravitational constant of main attractor (km^3 / s^2).
    r0 : Quantity
        Initial position (km).
    r : Quantity
        Final position (km).
    tof : Quantity
        Time of flight (s).
    M : int, optional
        Number of full revolutions, default to 0.
    numiter : int, optional
        Maximum number of iterations, default to 35.
    rtol : float, optional
        Relative tolerance of the algorithm, default to 1e-8.

    Yields
    ------
    v0, v : tuple
        Pair of velocity solutions.

    """
    k_ = k.to(u.km ** 3 / u.s ** 2).value
    r0_ = r0.to(u.km).value
    r_ = r.to(u.km).value
    tof_ = tof.to(u.s).value

    sols = _lambert_izzo(k_, r0_, r_, tof_, M, numiter, rtol)

    for v0, v in sols:
        yield v0 * u.km / u.s, v * u.km / u.s


def _lambert_izzo(k, r1, r2, tof, M, numiter, rtol):
    # Check preconditions
    assert tof > 0
    assert k > 0

    # Chord
    c = r2 - r1
    c_norm, r1_norm, r2_norm = norm(c), norm(r1), norm(r2)

    # Semiperimeter
    s = (r1_norm + r2_norm + c_norm) * .5

    # Versors
    i_r1, i_r2 = r1 / r1_norm, r2 / r2_norm
    i_h = np.cross(i_r1, i_r2)
    i_h = i_h / norm(i_h)  # Fixed from paper

    # Geometry of the problem
    ll = np.sqrt(1 - c_norm / s)

    if i_h[2] < 0:
        ll = -ll
        i_h = -i_h

    i_t1, i_t2 = np.cross(i_h, i_r1), np.cross(i_h, i_r2)  # Fixed from paper

    # Non dimensional time of flight
    T = np.sqrt(2 * k / s ** 3) * tof

    # Find solutions
    xy = _find_xy(ll, T, M, numiter, rtol)

    # Reconstruct
    gamma = np.sqrt(k * s / 2)
    rho = (r1_norm - r2_norm) / c_norm
    sigma = np.sqrt(1 - rho ** 2)

    for x, y in xy:
        v1, v2 = _reconstruct(x, y, r1_norm, r2_norm, i_r1, i_t1, i_r2, i_t2, ll, gamma, rho, sigma)
        yield v1, v2


def _reconstruct(x, y, r1, r2, i_r1, i_t1, i_r2, i_t2, ll, gamma, rho, sigma):
    """Reconstruct solution velocity vectors.

    """
    V_r1 = gamma * ((ll * y - x) - rho * (ll * y + x)) / r1
    V_r2 = -gamma * ((ll * y - x) + rho * (ll * y + x)) / r2
    V_t1 = gamma * sigma * (y + ll * x) / r1
    V_t2 = gamma * sigma * (y + ll * x) / r2
    v1 = V_r1 * i_r1 + V_t1 * i_t1
    v2 = V_r2 * i_r2 + V_t2 * i_t2
    return v1, v2


def _find_xy(ll, T, M, numiter, rtol):
    """Computes all x, y for given number of revolutions.

    """
    # For abs(ll) == 1 the derivative is not continuous
    assert abs(ll) < 1
    assert T > 0  # Mistake on original paper

    M_max = int(np.floor(T / pi))
    T_00 = np.arccos(ll) + ll * np.sqrt(1 - ll ** 2)  # T_xM

    # Refine maximum number of revolutions if necessary
    if T < T_00 + M_max * pi and M_max > 0:
        T_min = _compute_T_min(ll, M_max)
        if T < T_min:
            M_max -= 1

    # Check if a feasible solution exist for the given number of revolutions
    # This departs from the original paper in that we do not compute all solutions
    if M > M_max:
        raise ValueError("No feasible solution, try M <= {:d}".format(M_max))

    # Initial guess
    for x_0 in _initial_guess(T, ll, M):
        # Start Householder iterations from x_0 and find x, y
        x = householder(_tof_equation, x_0, args=(T, ll, M), tol=rtol, maxiter=numiter,
                        fprime=_tof_equation_p, fprime2=_tof_equation_pp,
                        fprime3=_tof_equation_ppp)
        y = _compute_y(x, ll)

        yield x, y


def _compute_y(x, ll):
    """Computes y.

    """
    return np.sqrt(1 - ll ** 2 * (1 - x ** 2))


def _compute_psi(x, ll):
    """Computes psi.

    "The auxiliary angle psi is computed using Eq.(17) by the appropriate
    inverse function"

    """
    y = _compute_y(x, ll)
    if -1 <= x < 1:
        # Elliptic motion
        # Use arc cosine to avoid numerical errors
        return np.arccos(x * y + ll * (1 - x ** 2))
    elif x > 1:
        # Hyperbolic motion
        # The hyperbolic sine is bijective
        return np.arcsinh((y - x * ll) * np.sqrt(x ** 2 - 1))
    else:
        # Parabolic motion
        return 0.0


def _tof_equation(x, T, ll, M):
    """Time of flight equation.

    """
    y = _compute_y(x, ll)
    if M == 0 and np.sqrt(0.6) < x < np.sqrt(1.4):
        eta = y - ll * x
        S_1 = (1 - ll - x * eta) * .5
        Q = 4 / 3 * hyp2f1(3, 1, 5 / 2, S_1)
        T_ = (eta ** 3 * Q + 4 * ll * eta) * .5
    else:
        psi = _compute_psi(x, ll)
        T_ = ((psi + M * pi) / np.sqrt(np.abs(1 - x ** 2)) - x + ll * y) / (1 - x ** 2)

    return T_ - T


def _tof_equation_p(x, _, ll, M):
    # TODO: What about derivatives when x approaches 1?
    y = _compute_y(x, ll)
    T = _tof_equation(x, 0.0, ll, M)
    return (3 * T * x - 2 + 2 * ll ** 3 * x / y) / (1 - x ** 2)


def _tof_equation_pp(x, _, ll, M):
    y = _compute_y(x, ll)
    T = _tof_equation(x, 0.0, ll, M)
    dT = _tof_equation_p(x, _, ll, M)
    return (3 * T + 5 * x * dT + 2 * (1 - ll ** 2) * ll ** 3 / y ** 3) / (1 - x ** 2)


def _tof_equation_ppp(x, _, ll, M):
    y = _compute_y(x, ll)
    T = _tof_equation(x, 0.0, ll, M)
    dT = _tof_equation_p(x, _, ll, M)
    ddT = _tof_equation_pp(x, _, ll, M)
    return (7 * x * ddT + 8 * dT - 6 * (1 - ll ** 2) * ll ** 5 * x / y ** 5) / (1 - x ** 2)


def _compute_T_min(ll, M):
    """Compute minimum T.

    """
    if ll == 1:
        x_T_min = 0.0
        T_min = _tof_equation(x_T_min, 0.0, ll, M)
    else:
        if M == 0:
            x_T_min = np.inf
            T_min = 0.0
        else:
            # Set x_0 > 0 to avoid problems at ll = -1
            x_T_min = newton(_tof_equation_p, 0.1, args=(0.0, ll, M),
                             fprime=_tof_equation_pp, fprime2=_tof_equation_ppp)
            T_min = _tof_equation(x_T_min, 0.0, ll, M)

    return x_T_min, T_min


def _initial_guess(T, ll, M):
    """Initial guess.

    """
    if M == 0:
        # Single revolution
        T_0 = np.arccos(ll) + ll * np.sqrt(1 - ll ** 2) + M * pi  # Equation 19
        T_1 = 2 * (1 - ll ** 3) / 3  # Equation 21
        if T >= T_0:
            x_0 = (T_0 / T) ** (2 / 3) - 1
        elif T < T_1:
            x_0 = 5 / 2 * T_1 / T * (T_1 - T) / (1 - ll ** 5) + 1
        else:
            # This is the real condition, which is not exactly equivalent
            # elif T_1 < T < T_0
            x_0 = (T_0 / T) ** (np.log2(T_1 / T_0)) - 1

        # Convert to tuple
        x_0 = x_0,
    else:
        # Multiple revolution
        x_0l = (((M * pi + pi) / (8 * T)) ** (2 / 3) - 1) / (((M * pi + pi) / (8 * T)) ** (2 / 3) + 1)
        x_0r = (((8 * T) / (M * pi)) ** (2 / 3) - 1) / (((8 * T) / (M * pi)) ** (2 / 3) + 1)

        x_0 = x_0l, x_0r

    return x_0


def householder(func, x0, fprime, fprime2, fprime3, args=(), tol=1e-8, maxiter=35):
    """Householder iteration scheme.

    """
    # FIXME: Use Halley's iteration for now
    return newton(func, x0, fprime, args, tol, maxiter, fprime2)
