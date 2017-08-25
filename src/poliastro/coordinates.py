"""Functions related to coordinate systems and transformations.

This module complements :py:mod:`astropy.coordinates`.

"""

import math

import astropy.time as time
import astropy.coordinates as coordinates
import astropy.units as u

from poliastro.util import transform
from poliastro.constants import J2000


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
        r = transform(r, -W, 'z')
        v = transform(v, -W, 'z')

    r_trans1 = transform(r, -(90 * u.deg - dec), 'x')
    r_trans2 = transform(r_trans1, -(90 * u.deg + ra), 'z')

    v_trans1 = transform(v, -(90 * u.deg - dec), 'x')
    v_trans2 = transform(v_trans1, -(90 * u.deg + ra), 'z')

    icrs_frame_pos_coord, icrs_frame_vel_coord = coordinates.get_body_barycentric_posvel(source_body.name, time=epoch)

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

    icrs_frame_pos_coord, icrs_frame_vel_coord = coordinates.get_body_barycentric_posvel(target_body.name, time=epoch)

    r_trans1 = r - icrs_frame_pos_coord.xyz
    r_trans2 = transform(r_trans1, (90 * u.deg + ra), 'z')
    r_f = transform(r_trans2, (90 * u.deg - dec), 'x')

    v_trans1 = v - icrs_frame_vel_coord.xyz
    v_trans2 = transform(v_trans1, (90 * u.deg + ra), 'z')
    v_f = transform(v_trans2, (90 * u.deg - dec), 'x')

    if rotate_meridian:
        r_f = transform(r_f, W, 'z')
        v_f = transform(v_f, W, 'z')

    return r_f.to(r.unit), v_f.to(v.unit)
