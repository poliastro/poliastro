"""Quasi optimal eccentricity-only change, with formulas developed by Pollard.

References
----------
* Pollard, J. E. "Simplified Approach for Assessment of Low-Thrust
  Elliptical Orbit Transfers", 1997.

"""
import numpy as np
from astropy import units as u
from numpy import cross
from numpy.linalg import norm

from poliastro.core.thrust.change_ecc_quasioptimal import extra_quantities


def change_ecc_quasioptimal(ss_0, ecc_f, f):
    """Guidance law from the model.
    Thrust is aligned with an inertially fixed direction perpendicular to the
    semimajor axis of the orbit.

    Parameters
    ----------
    ss_0 : Orbit
        Initial orbit, containing all the information.
    ecc_f : float
        Final eccentricity.
    f : float
        Magnitude of constant acceleration
    """
    # We fix the inertial direction at the beginning
    k = ss_0.attractor.k.to(u.km ** 3 / u.s ** 2).value
    a = ss_0.a.to(u.km).value
    ecc_0 = ss_0.ecc.value
    if ecc_0 > 0.001:  # Arbitrary tolerance
        ref_vec = ss_0.e_vec / ecc_0
    else:
        ref_vec = ss_0.r / norm(ss_0.r)

    h_unit = ss_0.h_vec / norm(ss_0.h_vec)
    thrust_unit = cross(h_unit, ref_vec) * np.sign(ecc_f - ecc_0)

    def a_d(t0, u_, k):
        accel_v = f * thrust_unit
        return accel_v

    delta_V, t_f = extra_quantities(k, a, ecc_0, ecc_f, f)
    return a_d, delta_V, t_f
