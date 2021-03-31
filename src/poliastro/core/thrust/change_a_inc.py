import numpy as np
from numba import njit as jit

from poliastro.core.util import circular_velocity


@jit
def extra_quantities(k, a_0, a_f, inc_0, inc_f, f):
    """Extra quantities given by the Edelbaum (a, i) model."""
    V_0, beta_0_, _ = compute_parameters(k, a_0, a_f, inc_0, inc_f)
    delta_V_ = delta_V(V_0, beta_0_, inc_0, inc_f)
    t_f_ = delta_V_ / f

    return delta_V_, t_f_


@jit
def beta(t, *, V_0, f, beta_0):
    """Compute yaw angle (β) as a function of time and the problem parameters."""
    return np.arctan2(V_0 * np.sin(beta_0), V_0 * np.cos(beta_0) - f * t)


@jit
def beta_0(V_0, V_f, inc_0, inc_f):
    """Compute initial yaw angle (β) as a function of the problem parameters."""
    delta_i_f = abs(inc_f - inc_0)
    return np.arctan2(
        np.sin(np.pi / 2 * delta_i_f), V_0 / V_f - np.cos(np.pi / 2 * delta_i_f)
    )


@jit
def compute_parameters(k, a_0, a_f, inc_0, inc_f):
    """Compute parameters of the model."""
    delta_inc = abs(inc_f - inc_0)
    V_0 = circular_velocity(k, a_0)
    V_f = circular_velocity(k, a_f)
    beta_0_ = beta_0(V_0, V_f, inc_0, inc_f)

    return V_0, beta_0_, delta_inc


@jit
def delta_V(V_0, beta_0, inc_0, inc_f):
    """Compute required increment of velocity."""
    delta_i_f = abs(inc_f - inc_0)
    return V_0 * np.cos(beta_0) - V_0 * np.sin(beta_0) / np.tan(
        np.pi / 2 * delta_i_f + beta_0
    )
