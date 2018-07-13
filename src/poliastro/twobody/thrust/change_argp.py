"""Argument of perigee change, with formulas developed by Pollard.
References
----------
* Pollard, J. E. "Simplified Approach for Assessment of Low-Thrust
  Elliptical Orbit Transfers", 1997.
* Pollard, J. E. "Evaluation of Low-Thrust Orbital Maneuvers", 1998.
"""

import numpy as np

from poliastro.twobody import rv
from poliastro.twobody.decorators import state_from_vector
from poliastro.util import norm_fast, circular_velocity_fast
from poliastro.jit import jit


@jit
def delta_V(V, ecc, argp_0, argp_f, f, A):
    """Compute required increment of velocity.
    """
    delta_argp = argp_f - argp_0
    return delta_argp / (3 * np.sign(delta_argp) / 2 * np.sqrt(1 - ecc ** 2) / ecc / V + A / f)


@jit
def extra_quantities(k, a, ecc, argp_0, argp_f, f, A=0.0):
    """Extra quantities given by the model.
    """
    V = circular_velocity_fast(k, a)
    delta_V_ = delta_V(V, ecc, argp_0, argp_f, f, A)
    t_f_ = delta_V_ / f

    return delta_V_, t_f_


def change_argp(k, a, ecc, argp_0, argp_f, f):
    """Guidance law from the model.
    Thrust is aligned with an inertially fixed direction perpendicular to the
    semimajor axis of the orbit.

    Parameters
    ----------
    f : float
        Magnitude of constant acceleration
    """

    def a_d(t0, u_, k):
        r = u_[:3]
        v = u_[3:]
        nu = rv.rv2coe(k, r, v)[-1]

        alpha_ = nu - np.pi / 2

        r_ = r / norm_fast(r)
        w_ = np.cross(r, v) / norm_fast(np.cross(r, v))
        s_ = np.cross(w_, r_)
        accel_v = f * (
            np.cos(alpha_) * s_ +
            np.sin(alpha_) * r_
        )
        return accel_v

    delta_V, t_f = extra_quantities(k, a, ecc, argp_0, argp_f, f, A=0.0)

    return a_d, delta_V, t_f
