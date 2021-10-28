import numpy as np
from numba import njit as jit
from numpy import cross
from numpy.linalg import norm

from poliastro.core.util import circular_velocity


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

    @jit
    def a_d(t0, u_, k):
        accel_v = f * thrust_unit
        return accel_v

    delta_V, t_f = extra_quantities(k, a, ecc_0, ecc_f, f)
    return a_d, delta_V, t_f
