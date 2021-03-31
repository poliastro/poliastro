import numpy as np
from numpy import cross
from numpy.linalg import norm

from poliastro.core.thrust.change_a_inc import (
    beta,
    compute_parameters,
    extra_quantities,
)


def change_a_inc(k, a_0, a_f, inc_0, inc_f, f):
    """Guidance law from the Edelbaum/Kéchichian theory, optimal transfer between circular inclined orbits
       (a_0, i_0) --> (a_f, i_f), ecc = 0.

    Parameters
    ----------
    k : float
        Gravitational parameter.
    a_0 : float
        Initial semimajor axis.
    a_f : float
        Final semimajor axis.
    inc_0 : float
        Initial inclination.
    inc_f : float
        Final inclination.
    f : float
        Magnitude of constant acceleration

    Notes
    -----
    Edelbaum theory, reformulated by Kéchichian.

    References
    ----------
    * Edelbaum, T. N. "Propulsion Requirements for Controllable
      Satellites", 1961.
    * Kéchichian, J. A. "Reformulation of Edelbaum's Low-Thrust
      Transfer Problem Using Optimal Control Theory", 1997.
    """

    V_0, beta_0_, _ = compute_parameters(k, a_0, a_f, inc_0, inc_f)

    def a_d(t0, u_, k):
        r = u_[:3]
        v = u_[3:]

        # Change sign of beta with the out-of-plane velocity
        beta_ = beta(t0, V_0=V_0, f=f, beta_0=beta_0_) * np.sign(r[0] * (inc_f - inc_0))

        t_ = v / norm(v)
        w_ = cross(r, v) / norm(cross(r, v))
        # n_ = cross(t_, w_)
        accel_v = f * (np.cos(beta_) * t_ + np.sin(beta_) * w_)
        return accel_v

    delta_V, t_f = extra_quantities(k, a_0, a_f, inc_0, inc_f, f)
    return a_d, delta_V, t_f
