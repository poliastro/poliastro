import numpy as np
from astropy import units as u

from poliastro.core.elements import (
    circular_velocity as circular_velocity_fast,
    coe2rv as coe2rv_fast,
    coe2rv_many as coe2rv_many_fast,
    eccentricity_vector as eccentricity_vector_fast,
)
from poliastro.core.propagation.farnocchia import (
    delta_t_from_nu as delta_t_from_nu_fast,
)

u_kms = u.km / u.s
u_km3s2 = u.km**3 / u.s**2


@u.quantity_input(k=u_km3s2, a=u.km)
def circular_velocity(k, a):
    """Circular velocity for a given body (k) and semimajor axis (a)."""
    return (
        circular_velocity_fast(k.to_value(u_km3s2), a.to_value(u.km)) * u_kms
    )


@u.quantity_input(k=u_km3s2, a=u.km)
def mean_motion(k, a):
    """Mean motion given body (k) and semimajor axis (a)."""
    return np.sqrt(k / abs(a**3)).to(1 / u.s) * u.rad


@u.quantity_input(k=u_km3s2, a=u.km)
def period(k, a):
    """Period given body (k) and semimajor axis (a)."""
    n = mean_motion(k, a)
    return 2 * np.pi * u.rad / n


@u.quantity_input(k=u_km3s2, r=u.km, v=u_kms)
def energy(k, r, v):
    """Specific energy."""
    return v @ v / 2 - k / np.sqrt(r @ r)


@u.quantity_input(k=u_km3s2, r=u.km, v=u_kms)
def eccentricity_vector(k, r, v):
    """Eccentricity vector."""
    return (
        eccentricity_vector_fast(
            k.to_value(u_km3s2), r.to_value(u.km), v.to_value(u_kms)
        )
        * u.one
    )


@u.quantity_input(nu=u.rad, ecc=u.one, k=u_km3s2, r_p=u.km)
def t_p(nu, ecc, k, r_p):
    """Elapsed time since latest perifocal passage."""
    # TODO: Make this a propagator method
    t_p = (
        delta_t_from_nu_fast(
            nu.to_value(u.rad),
            ecc.value,
            k.to_value(u_km3s2),
            r_p.to_value(u.km),
        )
        * u.s
    )
    return t_p


@u.quantity_input(
    k=u_km3s2,
    R=u.km,
    J2=u.one,
    n_sunsync=1 / u.s,
    a=u.km,
    ecc=u.one,
    inc=u.rad,
)
def heliosynchronous(k, R, J2, n_sunsync, a=None, ecc=None, inc=None):
    with np.errstate(invalid="raise"):
        if a is None and (ecc is not None) and (inc is not None):
            a = (
                -3
                * R**2
                * J2
                * np.sqrt(k)
                / (2 * n_sunsync * (1 - ecc**2) ** 2)
                * np.cos(inc)
            ) ** (2 / 7)
        elif ecc is None and (a is not None) and (inc is not None):
            ecc = np.sqrt(
                1
                - np.sqrt(
                    -3
                    * R**2
                    * J2
                    * np.sqrt(k)
                    * np.cos(inc)
                    / (2 * a ** (7 / 2) * n_sunsync)
                )
            )
        elif inc is None and (ecc is not None) and (a is not None):
            # Inclination is the unknown variable
            inc = np.arccos(
                -2
                * a ** (7 / 2)
                * n_sunsync
                * (1 - ecc**2) ** 2
                / (3 * R**2 * J2 * np.sqrt(k))
            )
        else:
            raise ValueError("Two parameters of (a, ecc, inc) are required")

    return a, ecc, inc


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
    ecc : ~astropy.units.Quantity, optional
        Eccentricity, default to None.

    """
    if ecc is None:
        ecc = 0.0549 * u.one

    return ecc


def coe2rv(k, p, ecc, inc, raan, argp, nu):
    rr, vv = coe2rv_fast(
        k.to_value(u_km3s2),
        p.to_value(u.km),
        ecc.to_value(u.one),
        inc.to_value(u.rad),
        raan.to_value(u.rad),
        argp.to_value(u.rad),
        nu.to_value(u.rad),
    )

    rr = rr << u.km
    vv = vv << (u.km / u.s)

    return rr, vv


def coe2rv_many(k_arr, p_arr, ecc_arr, inc_arr, raan_arr, argp_arr, nu_arr):
    rr_arr, vv_arr = coe2rv_many_fast(
        k_arr.to_value(u_km3s2),
        p_arr.to_value(u.km),
        ecc_arr.to_value(u.one),
        inc_arr.to_value(u.rad),
        raan_arr.to_value(u.rad),
        argp_arr.to_value(u.rad),
        nu_arr.to_value(u.rad),
    )

    rr_arr = rr_arr << u.km
    vv_arr = vv_arr << (u.km / u.s)

    return rr_arr, vv_arr
