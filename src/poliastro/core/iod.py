import numpy as np
from numba import njit as jit
from numpy import cross, pi
from numpy.linalg import norm

from poliastro.core.hyper import hyp2f1b
from poliastro.core.stumpff import c2, c3


@jit
def vallado(k, r0, r, tof, short, numiter, rtol):
    r"""Solves the Lambert's problem.

    The algorithm returns the initial velocity vector and the final one, these are
    computed by the following expresions:

    .. math::

        \vec{v_{o}} &= \frac{1}{g}(\vec{r} - f\vec{r_{0}}) \\
        \vec{v} &= \frac{1}{g}(\dot{g}\vec{r} - \vec{r_{0}})

    Therefore, the lagrange coefficients need to be computed. For the case of
    Lamber's problem, they can be expressed by terms of the initial and final vector:

    .. math::

        \begin{align}
            f = 1 -\frac{y}{r_{o}} \\
            g = A\sqrt{\frac{y}{\mu}} \\
            \dot{g} = 1 - \frac{y}{r} \\
        \end{align}

    Where y(z) is a function that depends on the :py:mod:`poliastro.core.stumpff` coefficients:

    .. math::

        y = r_{o} + r + A\frac{zS(z)-1}{\sqrt{C(z)}} \\
        A = \sin{(\Delta \nu)}\sqrt{\frac{rr_{o}}{1 - \cos{(\Delta \nu)}}}

    The value of z to evaluate the stump functions is solved by applying a Numerical method to
    the following equation:

    .. math::

        z_{i+1} = z_{i} - \frac{F(z_{i})}{F{}'(z_{i})}

    Function F(z)  to the expression:

    .. math::

        F(z) = \left [\frac{y(z)}{C(z)}  \right ]^{\frac{3}{2}}S(z) + A\sqrt{y(z)} - \sqrt{\mu}\Delta t

    Parameters
    ----------
    k: float
        Gravitational Parameter
    r0: ~np.array
        Initial position vector
    r: ~np.array
        Final position vector
    tof: ~float
        Time of flight
    numiter: int
        Number of iterations to
    rtol: int
        Number of revolutions

    Returns
    -------
    v0: ~np.array
        Initial velocity vector
    v: ~np.array
        Final velocity vector

    Examples
    --------

    >>> from poliastro.core.iod import vallado
    >>> from astropy import units as u
    >>> import numpy as np
    >>> from poliastro.bodies import Earth
    >>> k = Earth.k.to(u.km ** 3 / u.s ** 2)
    >>> r1 = np.array([5000, 10000, 2100]) * u.km # Initial position vector
    >>> r2 = np.array([-14600, 2500, 7000]) * u.km # Final position vector
    >>> tof = 3600 * u.s # Time of flight
    >>> v1, v2 = vallado(k.value, r1.value, r2.value, tof.value, short=True, numiter=35, rtol=1e-8)
    >>> v1 = v1 * u.km / u.s
    >>> v2 = v2 * u.km / u.s
    >>> print(v1, v2)
    [-5.99249503  1.92536671  3.24563805] km / s [-3.31245851 -4.19661901 -0.38528906] km / s

    Note
    ----
    This procedure can be found in section 5.3 of Curtis, with all the
    theoretical description of the problem. Analytical example can be found
    in the same book under name Example 5.2.

    """
    if short:
        t_m = +1
    else:
        t_m = -1

    norm_r0 = np.dot(r0, r0) ** 0.5
    norm_r = np.dot(r, r) ** 0.5
    norm_r0_times_norm_r = norm_r0 * norm_r
    norm_r0_plus_norm_r = norm_r0 + norm_r

    cos_dnu = np.dot(r0, r) / norm_r0_times_norm_r

    A = t_m * (norm_r * norm_r0 * (1 + cos_dnu)) ** 0.5

    if A == 0.0:
        raise RuntimeError("Cannot compute orbit, phase angle is 180 degrees")

    psi = 0.0
    psi_low = -4 * np.pi ** 2
    psi_up = 4 * np.pi ** 2

    count = 0

    while count < numiter:
        y = norm_r0_plus_norm_r + A * (psi * c3(psi) - 1) / c2(psi) ** 0.5
        if A > 0.0:
            # Readjust xi_low until y > 0.0
            # Translated directly from Vallado
            while y < 0.0:
                psi_low = psi
                psi = (
                    0.8
                    * (1.0 / c3(psi))
                    * (1.0 - norm_r0_times_norm_r * np.sqrt(c2(psi)) / A)
                )
                y = norm_r0_plus_norm_r + A * (psi * c3(psi) - 1) / c2(psi) ** 0.5

        xi = np.sqrt(y / c2(psi))
        tof_new = (xi ** 3 * c3(psi) + A * np.sqrt(y)) / np.sqrt(k)

        # Convergence check
        if np.abs((tof_new - tof) / tof) < rtol:
            break
        count += 1
        # Bisection check
        condition = tof_new <= tof
        psi_low = psi_low + (psi - psi_low) * condition
        psi_up = psi_up + (psi - psi_up) * (not condition)

        psi = (psi_up + psi_low) / 2
    else:
        raise RuntimeError("Maximum number of iterations reached")

    f = 1 - y / norm_r0
    g = A * np.sqrt(y / k)

    gdot = 1 - y / norm_r

    v0 = (r - f * r0) / g
    v = (gdot * r - r0) / g

    return v0, v


@jit
def izzo(k, r1, r2, tof, M, numiter, rtol):
    """Aplies izzo algorithm to solve Lambert's problem.

    Parameters
    ----------
    k: float
        Gravitational Constant
    r1: ~numpy.array
        Initial position vector
    r2: ~numpy.array
        Final position vector
    tof: float
        Time of flight between both positions
    M: int
        Number of revolutions
    numiter: int
        Numbert of iterations
    rtol: float
        Error tolerance

    Returns
    -------

    v1: ~numpy.array
        Initial velocity vector
    v2: ~numpy.array
        Final velocity vector

    """

    # Check preconditions
    assert tof > 0
    assert k > 0

    # Check collinearity of r1 and r2
    if np.all(cross(r1, r2) == 0):
        raise ValueError("Lambert solution cannot be computed for collinear vectors")

    # Chord
    c = r2 - r1
    c_norm, r1_norm, r2_norm = norm(c), norm(r1), norm(r2)

    # Semiperimeter
    s = (r1_norm + r2_norm + c_norm) * 0.5

    # Versors
    i_r1, i_r2 = r1 / r1_norm, r2 / r2_norm
    i_h = cross(i_r1, i_r2)
    i_h = i_h / norm(i_h)  # Fixed from paper

    # Geometry of the problem
    ll = np.sqrt(1 - min(1.0, c_norm / s))

    if i_h[2] < 0:
        ll = -ll
        i_h = -i_h

    i_t1, i_t2 = cross(i_h, i_r1), cross(i_h, i_r2)  # Fixed from paper

    # Non dimensional time of flight
    T = np.sqrt(2 * k / s ** 3) * tof

    # Find solutions
    xy = _find_xy(ll, T, M, numiter, rtol)

    # Reconstruct
    gamma = np.sqrt(k * s / 2)
    rho = (r1_norm - r2_norm) / c_norm
    sigma = np.sqrt(1 - rho ** 2)

    for x, y in xy:
        V_r1, V_r2, V_t1, V_t2 = _reconstruct(
            x, y, r1_norm, r2_norm, ll, gamma, rho, sigma
        )
        v1 = V_r1 * i_r1 + V_t1 * i_t1
        v2 = V_r2 * i_r2 + V_t2 * i_t2
        yield v1, v2


@jit
def _reconstruct(x, y, r1, r2, ll, gamma, rho, sigma):
    """Reconstruct solution velocity vectors."""
    V_r1 = gamma * ((ll * y - x) - rho * (ll * y + x)) / r1
    V_r2 = -gamma * ((ll * y - x) + rho * (ll * y + x)) / r2
    V_t1 = gamma * sigma * (y + ll * x) / r1
    V_t2 = gamma * sigma * (y + ll * x) / r2
    return [V_r1, V_r2, V_t1, V_t2]


@jit
def _find_xy(ll, T, M, numiter, rtol):
    """Computes all x, y for given number of revolutions."""
    # For abs(ll) == 1 the derivative is not continuous
    assert abs(ll) < 1
    assert T > 0  # Mistake on original paper

    M_max = np.floor(T / pi)
    T_00 = np.arccos(ll) + ll * np.sqrt(1 - ll ** 2)  # T_xM

    # Refine maximum number of revolutions if necessary
    if T < T_00 + M_max * pi and M_max > 0:
        _, T_min = _compute_T_min(ll, M_max, numiter, rtol)
        if T < T_min:
            M_max -= 1

    # Check if a feasible solution exist for the given number of revolutions
    # This departs from the original paper in that we do not compute all solutions
    if M > M_max:
        raise ValueError("No feasible solution, try lower M")

    # Initial guess
    for x_0 in _initial_guess(T, ll, M):
        # Start Householder iterations from x_0 and find x, y
        x = _householder(x_0, T, ll, M, rtol, numiter)
        y = _compute_y(x, ll)

        yield x, y


@jit
def _compute_y(x, ll):
    """Computes y."""
    return np.sqrt(1 - ll ** 2 * (1 - x ** 2))


@jit
def _compute_psi(x, y, ll):
    """Computes psi.

    "The auxiliary angle psi is computed using Eq.(17) by the appropriate
    inverse function"

    """
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


@jit
def _tof_equation(x, T0, ll, M):
    """Time of flight equation."""
    return _tof_equation_y(x, _compute_y(x, ll), T0, ll, M)


@jit
def _tof_equation_y(x, y, T0, ll, M):
    """Time of flight equation with externally computated y."""
    if M == 0 and np.sqrt(0.6) < x < np.sqrt(1.4):
        eta = y - ll * x
        S_1 = (1 - ll - x * eta) * 0.5
        Q = 4 / 3 * hyp2f1b(S_1)
        T_ = (eta ** 3 * Q + 4 * ll * eta) * 0.5
    else:
        psi = _compute_psi(x, y, ll)
        T_ = np.divide(
            np.divide(psi + M * pi, np.sqrt(np.abs(1 - x ** 2))) - x + ll * y,
            (1 - x ** 2),
        )

    return T_ - T0


@jit
def _tof_equation_p(x, y, T, ll):
    # TODO: What about derivatives when x approaches 1?
    return (3 * T * x - 2 + 2 * ll ** 3 * x / y) / (1 - x ** 2)


@jit
def _tof_equation_p2(x, y, T, dT, ll):
    return (3 * T + 5 * x * dT + 2 * (1 - ll ** 2) * ll ** 3 / y ** 3) / (1 - x ** 2)


@jit
def _tof_equation_p3(x, y, _, dT, ddT, ll):
    return (7 * x * ddT + 8 * dT - 6 * (1 - ll ** 2) * ll ** 5 * x / y ** 5) / (
        1 - x ** 2
    )


@jit
def _compute_T_min(ll, M, numiter, rtol):
    """Compute minimum T."""
    if ll == 1:
        x_T_min = 0.0
        T_min = _tof_equation(x_T_min, 0.0, ll, M)
    else:
        if M == 0:
            x_T_min = np.inf
            T_min = 0.0
        else:
            # Set x_i > 0 to avoid problems at ll = -1
            x_i = 0.1
            T_i = _tof_equation(x_i, 0.0, ll, M)
            x_T_min = _halley(x_i, T_i, ll, rtol, numiter)
            T_min = _tof_equation(x_T_min, 0.0, ll, M)

    return [x_T_min, T_min]


@jit
def _initial_guess(T, ll, M):
    """Initial guess."""
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

        return [x_0]
    else:
        # Multiple revolution
        x_0l = (((M * pi + pi) / (8 * T)) ** (2 / 3) - 1) / (
            ((M * pi + pi) / (8 * T)) ** (2 / 3) + 1
        )
        x_0r = (((8 * T) / (M * pi)) ** (2 / 3) - 1) / (
            ((8 * T) / (M * pi)) ** (2 / 3) + 1
        )

        return [x_0l, x_0r]


@jit
def _halley(p0, T0, ll, tol, maxiter):
    """Find a minimum of time of flight equation using the Halley method.

    Note
    ----
    This function is private because it assumes a calling convention specific to
    this module and is not really reusable.

    """
    for ii in range(maxiter):
        y = _compute_y(p0, ll)
        fder = _tof_equation_p(p0, y, T0, ll)
        fder2 = _tof_equation_p2(p0, y, T0, fder, ll)
        if fder2 == 0:
            raise RuntimeError("Derivative was zero")
        fder3 = _tof_equation_p3(p0, y, T0, fder, fder2, ll)

        # Halley step (cubic)
        p = p0 - 2 * fder * fder2 / (2 * fder2 ** 2 - fder * fder3)

        if abs(p - p0) < tol:
            return p
        p0 = p

    raise RuntimeError("Failed to converge")


@jit
def _householder(p0, T0, ll, M, tol, maxiter):
    """Find a zero of time of flight equation using the Householder method.

    Note
    ----
    This function is private because it assumes a calling convention specific to
    this module and is not really reusable.

    """
    for ii in range(maxiter):
        y = _compute_y(p0, ll)
        fval = _tof_equation_y(p0, y, T0, ll, M)
        T = fval + T0
        fder = _tof_equation_p(p0, y, T, ll)
        fder2 = _tof_equation_p2(p0, y, T, fder, ll)
        fder3 = _tof_equation_p3(p0, y, T, fder, fder2, ll)

        # Householder step (quartic)
        p = p0 - fval * (
            (fder ** 2 - fval * fder2 / 2)
            / (fder * (fder ** 2 - fval * fder2) + fder3 * fval ** 2 / 6)
        )

        if abs(p - p0) < tol:
            return p
        p0 = p

    raise RuntimeError("Failed to converge")
