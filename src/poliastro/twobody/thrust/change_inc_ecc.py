"""Simultaneous eccentricity and inclination changes.

References
----------
* Pollard, J. E. "Simplified Analysis of Low-Thrust Orbital Maneuvers", 2000.

"""
import numpy as np
from astropy import units as u
from numpy import cross
from numpy.linalg import norm

from poliastro.core.elements import rv2coe
from poliastro.core.thrust.change_inc_ecc import beta, extra_quantities


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

    # We fix the inertial direction at the beginning
    ecc_0 = ss_0.ecc.value
    if ecc_0 > 0.001:  # Arbitrary tolerance
        ref_vec = ss_0.e_vec / ecc_0
    else:
        ref_vec = ss_0.r / norm(ss_0.r)

    h_unit = ss_0.h_vec / norm(ss_0.h_vec)
    thrust_unit = cross(h_unit, ref_vec) * np.sign(ecc_f - ecc_0)

    inc_0 = ss_0.inc.to(u.rad).value
    argp = ss_0.argp.to(u.rad).value

    beta_0_ = beta(ecc_0, ecc_f, inc_0, inc_f, argp)

    def a_d(t0, u_, k):
        r = u_[:3]
        v = u_[3:]
        nu = rv2coe(k, r, v)[-1]
        beta_ = beta_0_ * np.sign(
            np.cos(nu)
        )  # The sign of ß reverses at minor axis crossings

        w_ = cross(r, v) / norm(cross(r, v))
        accel_v = f * (np.cos(beta_) * thrust_unit + np.sin(beta_) * w_)
        return accel_v

    delta_V, beta_, t_f = extra_quantities(
        ss_0.attractor.k.to(u.km ** 3 / u.s ** 2).value,
        ss_0.a.to(u.km).value,
        ecc_0,
        ecc_f,
        inc_0,
        inc_f,
        argp,
        f,
    )
    return a_d, delta_V, beta_, t_f
