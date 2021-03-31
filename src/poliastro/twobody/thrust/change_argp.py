"""Argument of perigee change, with formulas developed by Pollard.

References
----------
* Pollard, J. E. "Simplified Approach for Assessment of Low-Thrust
  Elliptical Orbit Transfers", 1997.
* Pollard, J. E. "Evaluation of Low-Thrust Orbital Maneuvers", 1998.

"""
import numpy as np
from numpy import cross
from numpy.linalg import norm

from poliastro.core.elements import rv2coe
from poliastro.core.thrust.change_argp import extra_quantities


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
        nu = rv2coe(k, r, v)[-1]

        alpha_ = nu - np.pi / 2

        r_ = r / norm(r)
        w_ = cross(r, v) / norm(cross(r, v))
        s_ = cross(w_, r_)
        accel_v = f * (np.cos(alpha_) * s_ + np.sin(alpha_) * r_)
        return accel_v

    delta_V, t_f = extra_quantities(k, a, ecc, argp_0, argp_f, f, A=0.0)

    return a_d, delta_V, t_f
