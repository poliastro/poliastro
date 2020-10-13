""" Low level propagation algorithms """

import numpy as np

from .danby import danby, danby_coe
from .farnocchia import farnocchia, farnocchia_coe
from .gooding import gooding, gooding_coe
from .markley import markley, markley_coe
from .mikkola import mikkola, mikkola_coe
from .pimienta import pimienta, pimienta_coe
from .vallado import vallado

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


def func_twobody(t0, u_, k, ad, ad_kwargs):
    """Differential equation for the initial value two body problem.

    This function follows Cowell's formulation.

    Parameters
    ----------
    t0 : float
        Time.
    u_ : ~numpy.ndarray
        Six component state vector [x, y, z, vx, vy, vz] (km, km/s).
    k : float
        Standard gravitational parameter.
    ad : function(t0, u, k)
        Non Keplerian acceleration (km/s2).
    ad_kwargs : optional
        perturbation parameters passed to ad.

    """
    ax, ay, az = ad(t0, u_, k, **ad_kwargs)

    x, y, z, vx, vy, vz = u_
    r3 = (x ** 2 + y ** 2 + z ** 2) ** 1.5

    du = np.array([vx, vy, vz, -k * x / r3 + ax, -k * y / r3 + ay, -k * z / r3 + az])
    return du
