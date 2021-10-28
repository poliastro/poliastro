"""Simultaneous eccentricity and inclination changes.

References
----------
* Pollard, J. E. "Simplified Analysis of Low-Thrust Orbital Maneuvers", 2000.

"""
from astropy import units as u

from poliastro.core.thrust.change_ecc_inc import change_ecc_inc as change_ecc_inc_fast


def change_ecc_inc(ss_0, ecc_f, inc_f, f):
    """Simultaneous eccentricity and inclination changes.
    Guidance law from the model.
    Thrust is aligned with an inertially fixed direction perpendicular to the
    semimajor axis of the orbit.

    Parameters
    ----------
    ss_0 : Orbit
        Initial orbit, containing all the information.
    ecc_f : ~astropy.units.quantity.Quantity
        Final eccentricity.
    inc_f : ~astropy.units.quantity.Quantity
        Final inclination (rad).
    f : ~astropy.units.quantity.Quantity
        Magnitude of constant acceleration (km / s**2).
    """

    a_d, delta_V, beta_, t_f = change_ecc_inc_fast(
        ss_0 = ss_0,
        ecc_f = ecc_f.to_value(u.one),
        inc_f = inc_f.to_value(u.rad),
        f = f.to_value(u.km / u.s ** 2),
    )

    return a_d, delta_V, beta_, t_f * u.s
