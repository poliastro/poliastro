import numpy as np
from numba import njit as jit

from poliastro._math.linalg import norm
from poliastro._math.special import stumpff_c2 as c2, stumpff_c3 as c3


@jit
def vallado(k, r0, v0, tof, numiter):
    r"""Solves Kepler's Equation by applying a Newton-Raphson method.

    If the position of a body along its orbit wants to be computed
    for a specific time, it can be solved by terms of the Kepler's Equation:

    .. math::
        E = M + e\sin{E}

    In this case, the equation is written in terms of the Universal Anomaly:

    .. math::

        \sqrt{\mu}\Delta t = \frac{r_{o}v_{o}}{\sqrt{\mu}}\chi^{2}C(\alpha \chi^{2}) + (1 - \alpha r_{o})\chi^{3}S(\alpha \chi^{2}) + r_{0}\chi

    This equation is solved for the universal anomaly by applying a Newton-Raphson numerical method.
    Once it is solved, the Lagrange coefficients are returned:

    .. math::

        \begin{align}
            f &= 1 \frac{\chi^{2}}{r_{o}}C(\alpha \chi^{2}) \\
            g &= \Delta t - \frac{1}{\sqrt{\mu}}\chi^{3}S(\alpha \chi^{2}) \\
            \dot{f} &= \frac{\sqrt{\mu}}{rr_{o}}(\alpha \chi^{3}S(\alpha \chi^{2}) - \chi) \\
            \dot{g} &= 1 - \frac{\chi^{2}}{r}C(\alpha \chi^{2}) \\
        \end{align}

    Lagrange coefficients can be related then with the position and velocity vectors:

    .. math::
        \begin{align}
            \vec{r} &= f\vec{r_{o}} + g\vec{v_{o}} \\
            \vec{v} &= \dot{f}\vec{r_{o}} + \dot{g}\vec{v_{o}} \\
        \end{align}

    Parameters
    ----------
    k : float
        Standard gravitational parameter.
    r0 : numpy.ndarray
        Initial position vector.
    v0 : numpy.ndarray
        Initial velocity vector.
    tof : float
        Time of flight.
    numiter : int
        Number of iterations.

    Returns
    -------
    f: float
        First Lagrange coefficient
    g: float
        Second Lagrange coefficient
    fdot: float
        Derivative of the first coefficient
    gdot: float
        Derivative of the second coefficient

    Notes
    -----
    The theoretical procedure is explained in section 3.7 of Curtis in really
    deep detail. For analytical example, check in the same book for example 3.6.

    """

    # Cache some results
    dot_r0v0 = r0 @ v0
    norm_r0 = norm(r0)
    sqrt_mu = k**0.5
    alpha = -(v0 @ v0) / k + 2 / norm_r0

    # First guess
    if alpha > 0:
        # Elliptic orbit
        xi_new = sqrt_mu * tof * alpha
    elif alpha < 0:
        # Hyperbolic orbit
        xi_new = (
            np.sign(tof)
            * (-1 / alpha) ** 0.5
            * np.log(
                (-2 * k * alpha * tof)
                / (
                    dot_r0v0
                    + np.sign(tof)
                    * np.sqrt(-k / alpha)
                    * (1 - norm_r0 * alpha)
                )
            )
        )
    else:
        # Parabolic orbit
        # (Conservative initial guess)
        xi_new = sqrt_mu * tof / norm_r0

    # Newton-Raphson iteration on the Kepler equation
    count = 0
    while count < numiter:
        xi = xi_new
        psi = xi * xi * alpha
        c2_psi = c2(psi)
        c3_psi = c3(psi)
        norm_r = (
            xi * xi * c2_psi
            + dot_r0v0 / sqrt_mu * xi * (1 - psi * c3_psi)
            + norm_r0 * (1 - psi * c2_psi)
        )
        xi_new = (
            xi
            + (
                sqrt_mu * tof
                - xi * xi * xi * c3_psi
                - dot_r0v0 / sqrt_mu * xi * xi * c2_psi
                - norm_r0 * xi * (1 - psi * c3_psi)
            )
            / norm_r
        )
        if abs(xi_new - xi) < 1e-7:
            break
        else:
            count += 1
    else:
        raise RuntimeError("Maximum number of iterations reached")

    # Compute Lagrange coefficients
    f = 1 - xi**2 / norm_r0 * c2_psi
    g = tof - xi**3 / sqrt_mu * c3_psi

    gdot = 1 - xi**2 / norm_r * c2_psi
    fdot = sqrt_mu / (norm_r * norm_r0) * xi * (psi * c3_psi - 1)

    return f, g, fdot, gdot
