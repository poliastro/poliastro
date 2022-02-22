import numpy as np
from numba import njit as jit
from numpy import cross

from poliastro._math.linalg import norm
from poliastro.core.elements import circular_velocity, rv2coe


@jit
def delta_V(V, ecc, argp_0, argp_f, f, A):
    """Compute required increment of velocity."""
    delta_argp = argp_f - argp_0
    return delta_argp / (
        3 * np.sign(delta_argp) / 2 * np.sqrt(1 - ecc**2) / ecc / V + A / f
    )


@jit
def extra_quantities(k, a, ecc, argp_0, argp_f, f, A=0.0):
    """Extra quantities given by the model."""
    V = circular_velocity(k, a)
    delta_V_ = delta_V(V, ecc, argp_0, argp_f, f, A)
    t_f_ = delta_V_ / f

    return delta_V_, t_f_


def change_argp(k, a, ecc, argp_0, argp_f, f):
    """Guidance law from the model.
    Thrust is aligned with an inertially fixed direction perpendicular to the
    semimajor axis of the orbit.

    Parameters
    ----------
    k : float
        Gravitational parameter (km**3 / s**2)
    a : float
        Semi-major axis (km)
    ecc : float
        Eccentricity
    argp_0 : float
        Initial argument of periapsis (rad)
    argp_f : float
        Final argument of periapsis (rad)
    f : float
        Magnitude of constant acceleration (km / s**2)

    Returns
    -------
    a_d : function
    delta_V : numpy.ndarray
    t_f : float
    """

    @jit
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
