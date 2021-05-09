import numpy as np
from astropy import units as u

from poliastro.core.util import circular_velocity as circular_velocity_fast

u_kms = u.km / u.s
u_km3s2 = u.km ** 3 / u.s ** 2


def circular_velocity(k, a):
    """Circular velocity for a given body (k) and semimajor axis (a)."""
    return circular_velocity_fast(k.to(u_km3s2).value, a.to(u.km).value) * u_kms


def mean_motion(k, a):
    """Mean motion given body (k) and semimajor axis (a)."""
    return np.sqrt(k / abs(a ** 3)).to(1 / u.s) * u.rad


def period(k, a):
    """Period given body (k) and semimajor axis (a)."""
    n = mean_motion(k, a)
    return 2 * np.pi * u.rad / n


@u.quantity_input(ecc=u.one)
def hyp_nu_limit(ecc, r_max_ratio=np.inf):
    r"""Limit true anomaly for hyperbolic orbits.

    Parameters
    ----------
    ecc : ~astropy.units.Quantity
        Eccentricity, should be larger than 1.
    r_max_ratio : float, optional
        Value of :math:`r_{\text{max}} / p` for this angle, default to infinity.

    """
    return np.arccos(-(1 - 1 / r_max_ratio) / ecc)


@u.quantity_input(R=u.m, J2=u.one, J3=u.one, a=u.m, inc=u.rad)
def get_eccentricity_critical_argp(R, J2, J3, a, inc):
    """Cccentricity for frozen orbits when the argument of perigee is critical.

    Parameters
    ----------
    R : ~astropy.units.Quantity
        Planet radius.
    J2 : ~astropy.units.Quantity
        Planet J2 coefficient.
    J3 : ~astropy.units.Quantity
        Planet J3 coefficient.
    a : ~astropy.units.Quantity
        Orbit's semimajor axis
    inc : ~astropy.units.Quantity
         Inclination.

    """
    ecc = -J3 * R * np.sin(inc) / 2 / J2 / a
    return ecc


@u.quantity_input(R=u.m, J2=u.one, J3=u.one, a=u.m, ecc=u.one)
def get_inclination_critical_argp(R, J2, J3, a, ecc):
    """Inclination for frozen orbits
    when the argument of perigee is critical and the eccentricity is given.

    Parameters
    ----------
    R : ~astropy.units.Quantity
        Planet radius.
    J2 : ~astropy.units.Quantity
        Planet J2 coefficient.
    J3 : ~astropy.units.Quantity
        Planet J3 coefficient.
    a : ~astropy.units.Quantity
        Semimajor axis.
    ecc : ~astropy.units.Quantity
        Eccentricity.

    """
    inc = np.arcsin(-ecc * a * J2 * 2 / R / J3) * u.rad
    return inc


@u.quantity_input(ecc=u.one)
def get_eccentricity_critical_inc(ecc=None):
    """Eccentricity for frozen orbits when the inclination is critical.

    If ecc is None we set an arbitrary value which is the Moon eccentricity
    because it seems reasonable.

    Parameters
    ----------
    ecc: ~astropy.units.Quantity, optional
        Eccentricity, default to None.

    """
    if ecc is None:
        ecc = 0.0549 * u.one

    return ecc
