import numpy as np
from numba import njit as jit

from poliastro.core.util import circular_velocity


@jit
def delta_V(V, ecc, argp_0, argp_f, f, A):
    """Compute required increment of velocity."""
    delta_argp = argp_f - argp_0
    return delta_argp / (
        3 * np.sign(delta_argp) / 2 * np.sqrt(1 - ecc ** 2) / ecc / V + A / f
    )


@jit
def extra_quantities(k, a, ecc, argp_0, argp_f, f, A=0.0):
    """Extra quantities given by the model."""
    V = circular_velocity(k, a)
    delta_V_ = delta_V(V, ecc, argp_0, argp_f, f, A)
    t_f_ = delta_V_ / f

    return delta_V_, t_f_
