"""Quasi optimal eccentricity-only change, with formulas developed by Pollard.

References
----------
* Pollard, J. E. "Simplified Approach for Assessment of Low-Thrust
  Elliptical Orbit Transfers", 1997.

"""
from astropy import units as u

from poliastro.core.thrust.change_ecc_quasioptimal import change_ecc_quasioptimal as change_ecc_quasioptimal_fast


def change_ecc_quasioptimal(ss_0, ecc_f, f):
    """Guidance law from the model.
    Thrust is aligned with an inertially fixed direction perpendicular to the
    semimajor axis of the orbit.

    Parameters
    ----------
    ss_0 : Orbit
        Initial orbit, containing all the information.
    ecc_f : ~astropy.units.quantity.Quantity
        Final eccentricity.
    f : ~astropy.units.quantity.Quantity
        Magnitude of constant acceleration (km / s**2).
    """
    a_d, delta_V, t_f = change_ecc_quasioptimal_fast(
        ss_0 = ss_0,
        ecc_f = ecc_f.to_value(u.one),
        f = f.to_value(u.km / u.s ** 2),
    )

    return a_d, delta_V, t_f * u.s
