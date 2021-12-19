import numpy as np
from numba import njit as jit

from poliastro.core.elements import circular_velocity


@jit
def delta_V(V_0, ecc_0, ecc_f):
    """Compute required increment of velocity."""
    return 2 / 3 * V_0 * np.abs(np.arcsin(ecc_0) - np.arcsin(ecc_f))


@jit
def extra_quantities(k, a, ecc_0, ecc_f, f):
    """Extra quantities given by the model."""
    V_0 = circular_velocity(k, a)
    delta_V_ = delta_V(V_0, ecc_0, ecc_f)
    t_f_ = delta_V_ / f

    return delta_V_, t_f_
