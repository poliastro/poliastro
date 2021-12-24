"""Low level propagation algorithms."""

import numpy as np
from numba import njit as jit

from poliastro.core.propagation.danby import danby, danby_coe
from poliastro.core.propagation.farnocchia import farnocchia, farnocchia_coe
from poliastro.core.propagation.gooding import gooding, gooding_coe
from poliastro.core.propagation.markley import markley, markley_coe
from poliastro.core.propagation.mikkola import mikkola, mikkola_coe
from poliastro.core.propagation.pimienta import pimienta, pimienta_coe
from poliastro.core.propagation.vallado import vallado

__all__ = [
    "func_twobody",
    "farnocchia_coe",
    "farnocchia",
    "vallado",
    "mikkola_coe",
    "mikkola",
    "markley_coe",
    "markley",
    "pimienta_coe",
    "pimienta",
    "gooding_coe",
    "gooding",
    "danby_coe",
    "danby",
]


@jit
def func_twobody(t0, u_, k):
    """Differential equation for the initial value two body problem.

    This function follows Cowell's formulation.

    Parameters
    ----------
    t0 : float
        Time.
    u_ : ~np.array
        Six component state vector [x, y, z, vx, vy, vz] (km, km/s).
    k : float
        Standard gravitational parameter.

    """
    x, y, z, vx, vy, vz = u_
    r3 = (x ** 2 + y ** 2 + z ** 2) ** 1.5

    du = np.array([vx, vy, vz, -k * x / r3, -k * y / r3, -k * z / r3])
    return du
