"""Simultaneous eccentricity and inclination changes.
References
----------
* Pollard, J. E. "Simplified Analysis of Low-Thrust Orbital Maneuvers", 2000.
"""
import numpy as np

from astropy import units as u

from poliastro.twobody.decorators import state_from_vector
from poliastro.util import norm, circular_velocity
from poliastro.twobody import rv
from poliastro.jit import jit


@jit
def beta(ecc_0, ecc_f, inc_0, inc_f, argp):
    # Note: "The argument of perigee will vary during the orbit transfer
    # due to the natural drift and because e may approach zero.
    # However, [the equation] still gives a good estimate of the desired
    # thrust angle."
    return np.arctan(abs(3 * np.pi * (inc_f - inc_0) / (4 * np.cos(argp) * (ecc_0 - ecc_f + np.log(
        (1 + ecc_f) * (-1 + ecc_0) / ((1 + ecc_0) * (-1 + ecc_f)))))))


@jit
def delta_V(V_0, ecc_0, ecc_f, beta_):
    """Compute required increment of velocity.
    """
    return 2 * V_0 * np.abs(np.arcsin(ecc_0) - np.arcsin(ecc_f)) / (3 * np.cos(beta_))


def extra_quantities(k, a, ecc_0, ecc_f, inc_0, inc_f, argp, f):
    """Extra quantities given by the model.
    """
    beta_ = beta(ecc_0, ecc_f, inc_0, inc_f, argp)
    V_0 = circular_velocity(k, a)
    delta_V_ = delta_V(V_0, ecc_0, ecc_f, beta_)
    t_f_ = delta_V_ / f

    return delta_V_, beta_, t_f_


def change_inc_ecc(ss_0, ecc_f, inc_f, f):
    """Guidance law from the model.
    Thrust is aligned with an inertially fixed direction perpendicular to the
    semimajor axis of the orbit.
    Parameters
    ----------
    ss_0 : Orbit
        Initial orbit, containing all the information.
    ecc_f : float
        Final eccentricity.
    inc_f : float
        Final inclination.
    f : float
        Magnitude of constant acceleration.
    """
    k = ss_0.attractor.k.to(u.km**3 / u.s**2).value
    a = ss_0.a.to(u.km).value
    # We fix the inertial direction at the beginning
    ecc_0 = ss_0.ecc.value
    if ecc_0 > 0.001:  # Arbitrary tolerance
        ref_vec = ss_0.e_vec / ecc_0
    else:
        ref_vec = ss_0.r / norm(ss_0.r)

    h_unit = ss_0.h_vec / norm(ss_0.h_vec)
    thrust_unit = np.cross(h_unit, ref_vec) * np.sign(ecc_f - ecc_0)

    inc_0 = ss_0.inc.to(u.rad).value
    argp = ss_0.argp.to(u.rad).value

    beta_0_ = beta(ecc_0, ecc_f, inc_0, inc_f, argp)

    def a_d(t0, u_, k):
        r = u_[:3]
        v = u_[3:]
        nu = rv.rv2coe(k, r, v)[-1]
        beta_ = beta_0_ * np.sign(np.cos(nu))  # The sign of ß reverses at minor axis crossings

        w_ = np.cross(r, v) / norm(np.cross(r, v))
        accel_v = f * (
            np.cos(beta_) * thrust_unit +
            np.sin(beta_) * w_
        )
        return accel_v

    delta_V, beta_, t_f = extra_quantities(k, a, ecc_0, ecc_f, inc_0, inc_f, argp, f)
    return a_d, delta_V, beta_, t_f
