"""Simultaneous eccentricity and inclination changes.

References
----------
* Pollard, J. E. "Simplified Analysis of Low-Thrust Orbital Maneuvers", 2000.

"""
import numpy as np
from numba import njit as jit
from numpy import cross

from poliastro._math.linalg import norm
from poliastro.core.elements import (
    circular_velocity,
    eccentricity_vector,
    rv2coe,
)


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
                    + np.log(
                        (1 + ecc_f)
                        * (-1 + ecc_0)
                        / ((1 + ecc_0) * (-1 + ecc_f))
                    )
                )
            )
        )
    )


@jit
def delta_V(V_0, ecc_0, ecc_f, beta_):
    """Compute required increment of velocity."""
    return (
        2
        * V_0
        * np.abs(np.arcsin(ecc_0) - np.arcsin(ecc_f))
        / (3 * np.cos(beta_))
    )


@jit
def delta_t(delta_v, f):
    """Compute required increment of velocity."""
    return delta_v / f


def change_ecc_inc(k, a, ecc_0, ecc_f, inc_0, inc_f, argp, r, v, f):
    # We fix the inertial direction at the beginning
    if ecc_0 > 0.001:  # Arbitrary tolerance
        e_vec = eccentricity_vector(k, r, v)
        ref_vec = e_vec / ecc_0
    else:
        ref_vec = r / norm(r)

    h_vec = cross(r, v)  # Specific angular momentum vector
    h_unit = h_vec / norm(h_vec)
    thrust_unit = cross(h_unit, ref_vec) * np.sign(ecc_f - ecc_0)

    beta_0 = beta(ecc_0, ecc_f, inc_0, inc_f, argp)

    @jit
    def a_d(t0, u_, k_):
        r_ = u_[:3]
        v_ = u_[3:]
        nu = rv2coe(k_, r_, v_)[-1]
        beta_ = beta_0 * np.sign(
            np.cos(nu)
        )  # The sign of ß reverses at minor axis crossings

        w_ = cross(r_, v_) / norm(cross(r_, v_))
        accel_v = f * (np.cos(beta_) * thrust_unit + np.sin(beta_) * w_)
        return accel_v

    delta_v = delta_V(circular_velocity(k, a), ecc_0, ecc_f, beta_0)
    t_f = delta_t(delta_v, f)

    return a_d, delta_v, t_f
