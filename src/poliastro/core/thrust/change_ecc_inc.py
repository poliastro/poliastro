import numpy as np
from numba import njit as jit

from poliastro.core.util import circular_velocity


@jit
def beta(ecc_0, ecc_f, inc_0, inc_f, argp):
    # Note: "The argument of perigee will vary during the orbit transfer
    # due to the natural drift and because e may approach zero.
    # However, [the equation] still gives a good estimate of the desired
    # thrust angle."
    return np.arctan(
        abs(
            3
            * np.pi
            * (inc_f - inc_0)
            / (
                4
                * np.cos(argp)
                * (
                    ecc_0
                    - ecc_f
                    + np.log((1 + ecc_f) * (-1 + ecc_0) / ((1 + ecc_0) * (-1 + ecc_f)))
                )
            )
        )
    )


@jit
def delta_V(V_0, ecc_0, ecc_f, beta_):
    """Compute required increment of velocity."""
    return 2 * V_0 * np.abs(np.arcsin(ecc_0) - np.arcsin(ecc_f)) / (3 * np.cos(beta_))


@jit
def extra_quantities(k, a, ecc_0, ecc_f, inc_0, inc_f, argp, f):
    """Extra quantities given by the model."""
    beta_ = beta(ecc_0, ecc_f, inc_0, inc_f, argp)
    V_0 = circular_velocity(k, a)
    delta_V_ = delta_V(V_0, ecc_0, ecc_f, beta_)
    t_f_ = delta_V_ / f

    return delta_V_, beta_, t_f_
