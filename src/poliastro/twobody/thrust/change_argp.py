"""Argument of perigee change, with formulas developed by Pollard.

References
----------
* Pollard, J. E. "Simplified Approach for Assessment of Low-Thrust
  Elliptical Orbit Transfers", 1997.
* Pollard, J. E. "Evaluation of Low-Thrust Orbital Maneuvers", 1998.

"""
from astropy import units as u

from poliastro.core.thrust.change_argp import change_argp as change_a_inc_fast


def change_argp(k, a, ecc, argp_0, argp_f, f):
    """Guidance law from the model.
    Thrust is aligned with an inertially fixed direction perpendicular to the
    semimajor axis of the orbit.

    Parameters
    ----------
    f : ~astropy.units.quantity.Quantity
        Magnitude of constant acceleration
    """
    a_d, delta_V, t_f = change_a_inc_fast(
        k,
        a,
        ecc,
        argp_0,
        argp_f,
        f=f.to_value(u.km / u.s ** 2),
    )

    return a_d, delta_V, t_f * u.s
