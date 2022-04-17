import numpy as np
from numba import njit as jit
from numpy import cross

from poliastro._math.linalg import norm
from poliastro.core.elements import circular_velocity


@jit
def extra_quantities(k, a_0, a_f, inc_0, inc_f, f):
    """Extra quantities given by the Edelbaum (a, i) model."""
    V_0, V_f, beta_0_ = compute_parameters(k, a_0, a_f, inc_0, inc_f)
    delta_V_ = delta_V(V_0, V_f, beta_0_, inc_0, inc_f)
    t_f_ = delta_V_ / f

    return delta_V_, t_f_


@jit
def beta(t, V_0, f, beta_0):
    """Compute yaw angle (β) as a function of time and the problem parameters."""
    return np.arctan2(V_0 * np.sin(beta_0), V_0 * np.cos(beta_0) - f * t)


@jit
def beta_0(V_0, V_f, inc_0, inc_f):
    """Compute initial yaw angle (β) as a function of the problem parameters."""
    delta_i_f = abs(inc_f - inc_0)
    return np.arctan2(
        np.sin(np.pi / 2 * delta_i_f),
        V_0 / V_f - np.cos(np.pi / 2 * delta_i_f),
    )


@jit
def compute_parameters(k, a_0, a_f, inc_0, inc_f):
    """Compute parameters of the model."""
    V_0 = circular_velocity(k, a_0)
    V_f = circular_velocity(k, a_f)
    beta_0_ = beta_0(V_0, V_f, inc_0, inc_f)

    return V_0, V_f, beta_0_


@jit
def delta_V(V_0, V_f, beta_0, inc_0, inc_f):
    """Compute required increment of velocity."""
    delta_i_f = abs(inc_f - inc_0)
    if delta_i_f == 0:
        return abs(V_f - V_0)
    return V_0 * np.cos(beta_0) - V_0 * np.sin(beta_0) / np.tan(
        np.pi / 2 * delta_i_f + beta_0
    )


def change_a_inc(k, a_0, a_f, inc_0, inc_f, f):
    """Change semimajor axis and inclination.
       Guidance law from the Edelbaum/Kéchichian theory, optimal transfer between circular inclined orbits
       (a_0, i_0) --> (a_f, i_f), ecc = 0.

    Parameters
    ----------
    k : float
        Gravitational parameter.
    a_0 : float
        Initial semimajor axis (km).
    a_f : float
        Final semimajor axis (km).
    inc_0 : float
        Initial inclination (rad).
    inc_f : float
        Final inclination (rad).
    f : float
        Magnitude of constant acceleration (km / s**2).

    Returns
    -------
    a_d : function
    delta_V : numpy.ndarray
    t_f : float

    Notes
    -----
    Edelbaum theory, reformulated by Kéchichian.

    References
    ----------
    * Edelbaum, T. N. "Propulsion Requirements for Controllable
      Satellites", 1961.
    * Kéchichian, J. A. "Reformulation of Edelbaum's Low-Thrust
      Transfer Problem Using Optimal Control Theory", 1997.
    """

    V_0, V_f, beta_0_ = compute_parameters(k, a_0, a_f, inc_0, inc_f)

    @jit
    def a_d(t0, u_, k):
        r = u_[:3]
        v = u_[3:]

        # Change sign of beta with the out-of-plane velocity
        beta_ = beta(t0, V_0, f, beta_0_) * np.sign(r[0] * (inc_f - inc_0))

        t_ = v / norm(v)
        w_ = cross(r, v) / norm(cross(r, v))
        accel_v = f * (np.cos(beta_) * t_ + np.sin(beta_) * w_)
        return accel_v

    delta_V, t_f = extra_quantities(k, a_0, a_f, inc_0, inc_f, f)
    return a_d, delta_V, t_f
