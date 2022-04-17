from astropy import units as u

from poliastro.core.thrust.change_a_inc import (
    change_a_inc as change_a_inc_fast,
)


def change_a_inc(k, a_0, a_f, inc_0, inc_f, f):
    """Change semimajor axis and inclination.
       Guidance law from the Edelbaum/Kéchichian theory, optimal transfer between circular inclined orbits
       (a_0, i_0) --> (a_f, i_f), ecc = 0.

    Parameters
    ----------
    k : ~astropy.units.quantity.Quantity
        Gravitational parameter.
    a_0 : ~astropy.units.quantity.Quantity
        Initial semimajor axis (km).
    a_f : ~astropy.units.quantity.Quantity
        Final semimajor axis (km).
    inc_0 : ~astropy.units.quantity.Quantity
        Initial inclination (rad).
    inc_f : ~astropy.units.quantity.Quantity
        Final inclination (rad).
    f : ~astropy.units.quantity.Quantity
        Magnitude of constant acceleration (km / s**2).

    Returns
    -------
    a_d, delta_V, t_f : tuple (function, ~astropy.units.quantity.Quantity, ~astropy.time.Time)

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
    a_d, delta_V, t_f = change_a_inc_fast(
        k=k.to_value(u.km**3 / u.s**2),
        a_0=a_0.to_value(u.km),
        a_f=a_f.to_value(u.km),
        inc_0=inc_0.to_value(u.rad),
        inc_f=inc_f.to_value(u.rad),
        f=f.to_value(u.km / u.s**2),
    )
    return a_d, delta_V, t_f * u.s
