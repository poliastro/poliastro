"""Functions related to coordinate systems and transformations.

This module complements :py:mod:`astropy.coordinates`.

"""
from math import cos, sin, sqrt

import astropy.units as u
import numpy as np
from astropy.coordinates import get_body_barycentric_posvel

from poliastro.constants import J2000
from poliastro.core.elements import rv2coe
from poliastro.util import transform as transform_vector


def body_centered_to_icrs(r, v, source_body, epoch=J2000, rotate_meridian=False):
    """Converts position and velocity body-centered frame to ICRS.

    Parameters
    ----------
    r : ~astropy.units.Quantity
        Position vector in a body-centered reference frame.
    v : ~astropy.units.Quantity
        Velocity vector in a body-centered reference frame.
    source_body : Body
        Source body.
    epoch : ~astropy.time.Time, optional
        Epoch, default to J2000.
    rotate_meridian : bool, optional
        Whether to apply the rotation of the meridian too, default to False.

    Returns
    -------
    r, v : tuple (~astropy.units.Quantity)
        Position and velocity vectors in ICRS.


    """

    ra, dec, W = source_body.rot_elements_at_epoch(epoch)
    if rotate_meridian:
        r = transform_vector(r, -W, "z")
        v = transform_vector(v, -W, "z")

    r_trans1 = transform_vector(r, -(90 * u.deg - dec), "x")
    r_trans2 = transform_vector(r_trans1, -(90 * u.deg + ra), "z")

    v_trans1 = transform_vector(v, -(90 * u.deg - dec), "x")
    v_trans2 = transform_vector(v_trans1, -(90 * u.deg + ra), "z")

    icrs_frame_pos_coord, icrs_frame_vel_coord = get_body_barycentric_posvel(
        source_body.name, time=epoch
    )

    r_f = icrs_frame_pos_coord.xyz + r_trans2
    v_f = icrs_frame_vel_coord.xyz + v_trans2

    return r_f.to(r.unit), v_f.to(v.unit)


def icrs_to_body_centered(r, v, target_body, epoch=J2000, rotate_meridian=False):
    """Converts position and velocity in ICRS to body-centered frame.

    Parameters
    ----------
    r : ~astropy.units.Quantity
        Position vector in ICRS.
    v : ~astropy.units.Quantity
        Velocity vector in ICRS.
    target_body : Body
        Target body.
    epoch : ~astropy.time.Time, optional
        Epoch, default to J2000.
    rotate_meridian : bool, optional
        Whether to apply the rotation of the meridian too, default to False.

    Returns
    -------
    r, v : tuple (~astropy.units.Quantity)
        Position and velocity vectors in a body-centered reference frame.

    """

    ra, dec, W = target_body.rot_elements_at_epoch(epoch)

    icrs_frame_pos_coord, icrs_frame_vel_coord = get_body_barycentric_posvel(
        target_body.name, time=epoch
    )

    r_trans1 = r - icrs_frame_pos_coord.xyz
    r_trans2 = transform_vector(r_trans1, (90 * u.deg + ra), "z")
    r_f = transform_vector(r_trans2, (90 * u.deg - dec), "x")

    v_trans1 = v - icrs_frame_vel_coord.xyz
    v_trans2 = transform_vector(v_trans1, (90 * u.deg + ra), "z")
    v_f = transform_vector(v_trans2, (90 * u.deg - dec), "x")

    if rotate_meridian:
        r_f = transform_vector(r_f, W, "z")
        v_f = transform_vector(v_f, W, "z")

    return r_f.to(r.unit), v_f.to(v.unit)


def inertial_body_centered_to_pqw(r, v, source_body):
    """Converts position and velocity from inertial body-centered frame to perifocal frame.

    Parameters
    ----------
    r : ~astropy.units.Quantity
        Position vector in a inertial body-centered reference frame.
    v : ~astropy.units.Quantity
        Velocity vector in a inertial body-centered reference frame.
    source_body : Body
        Source body.

    Returns
    -------
    r_pqw, v_pqw : tuple (~astropy.units.Quantity)
        Position and velocity vectors in ICRS.


    """
    r = r.to("km").value
    v = v.to("km/s").value
    k = source_body.k.to("km^3 / s^2").value

    p, ecc, inc, _, _, nu = rv2coe(k, r, v)

    r_pqw = (np.array([cos(nu), sin(nu), 0 * nu]) * p / (1 + ecc * cos(nu))).T * u.km
    v_pqw = (np.array([-sin(nu), (ecc + cos(nu)), 0]) * sqrt(k / p)).T * u.km / u.s

    return r_pqw, v_pqw
