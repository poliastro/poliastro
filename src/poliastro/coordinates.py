"""Functions related to coordinate systems and transformations.

This module complements :py:mod:`astropy.coordinates`.

"""

import math

import astropy.time as time
import astropy.coordinates as coordinates
import astropy.units as u

from poliastro.util import transform


T = 0
d = 0
J2000 = time.Time('J2000', scale='tdb')

Ja = 99.360714 + 4850.4046 * T
Jb = 175.895369 + 1191.9605 * T
Jc = 300.323162 + 262.5475 * T
Jd = 114.012305 + 6070.2476 * T
Je = 49.511251 + 64.3000 * T

N = 357.85 + 52.316 * T

bodies_orientations = {
    "Mercury": {
        "ra": 281.01 - 0.033 * T,
        "dec": 61.45 - 0.005 * T,
        "W": 329.548 + 6.1385025 * d
    },
    "Venus": {
        "ra": 272.76,
        "dec": 67.16,
        "W": 160.20 - 1.4813688 * d
    },
    "Mars": {
        "ra": 317.68143 - 0.1061 * T,
        "dec": 52.88650 - 0.0609 * T,
        "W": 176.630 + 350.89198226 * d
    },
    "Jupiter": {
        "ra": 268.056595 - 0.006499 * T + 0.000117 * math.sin(Ja) + 0.000938 * math.sin(Jb)
              + 0.001432 * math.sin(Jc) + 0.000030 * math.sin(Jd) + 0.002150 * math.sin(Je),
        "dec": 64.495303 + 0.002413 * T + 0.000050 * math.cos(Ja) + 0.000404 * math.cos(Jb)
               + 0.000617 * math.cos(Jc) - 0.000013 * math.cos(Jd) + 0.000926 * math.cos(Je),
        "W": 284.95 + 870.5366420 * d
    },
    "Saturn": {
        "ra": 40.589 - 0.036 * T,
        "dec": 83.537 - 0.004 * T,
        "W": 38.90 + 810.7939024 * d
    },
    "Uranus": {
        "ra": 257.311,
        "dec": -15.175,
        "W": 203.81 - 501.1600928 * d
    },
    "Neptune": {
        "ra": 299.36 + 0.70 * math.sin(N),
        "dec": 43.46 - 0.51 * math.cos(N),
        "W": 253.18 + 536.3128492 * d - 0.48 * math.sin(N)
    },
    "Pluto": {
        "ra": 312.993,
        "dec": 6.163,
        "W": 237.305 - 56.3625225 * d
    },
    "Earth": {
        "ra": 0.00 - 0.641 * T,
        "dec": 90.00 - 0.557 * T,
        "W": 190.147 + 360.9856235 * d
    }
}


def body_centered_to_icrs(r, v, source_body, epoch=None, rotate_meridian=False):
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

    # All epochs must be in TDB
    if not epoch:
        epoch = time.Time('J2000', scale='tdb')

    T = (epoch.tdb - J2000).to('day').value / 36525
    d = (epoch.tdb - J2000).to('day').value

    if rotate_meridian:
        r = transform(r, -bodies_orientations[source_body.name]['W'], 'z', u.deg)
        v = transform(v, -bodies_orientations[source_body.name]['W'], 'z', u.deg)

    r_trans1 = transform(r, -(90 - bodies_orientations[source_body.name]['dec']), 'x', u.deg)
    r_trans2 = transform(r_trans1, -(90 + bodies_orientations[source_body.name]['ra']), 'z', u.deg)

    v_trans1 = transform(v, -(90 - bodies_orientations[source_body.name]['dec']), 'x', u.deg)
    v_trans2 = transform(v_trans1, -(90 + bodies_orientations[source_body.name]['ra']), 'z', u.deg)

    icrs_frame_pos_coord, icrs_frame_vel_coord = coordinates.get_body_barycentric_posvel(source_body.name, time=epoch)

    r_f = icrs_frame_pos_coord.xyz + r_trans2
    v_f = icrs_frame_vel_coord.xyz + v_trans2

    # icrs_coord = coordinates.ICRS(x=r_f[0], y=r_f[1], z=r_f[2], v_x=v_f[0], v_y=v_f[1],
    #                              v_z=v_f[2], representation='cartesian')
    return r_f, v_f


def icrs_to_body_centered(r, v, target_body, epoch=None, rotate_meridian=False):
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

    T = (epoch.tdb - J2000).to('day').value / 36525
    d = (epoch.tdb - J2000).to('day').value

    icrs_frame_pos_coord, icrs_frame_vel_coord = coordinates.get_body_barycentric_posvel(target_body.name, time=epoch)

    r_trans1 = r - icrs_frame_pos_coord.xyz
    r_trans2 = transform(r_trans1, (90 + bodies_orientations[target_body.name]['ra']), 'z', u.deg)
    r_f = transform(r_trans2, (90 - bodies_orientations[target_body.name]['dec']), 'x', u.deg)

    v_trans1 = v - icrs_frame_vel_coord.xyz
    v_trans2 = transform(v_trans1, (90 + bodies_orientations[target_body.name]['ra']), 'z', u.deg)
    v_f = transform(v_trans2, (90 - bodies_orientations[target_body.name]['dec']), 'x', u.deg)

    if rotate_meridian:
        r_f = transform(r_f, bodies_orientations[target_body.name]['W'], 'z', u.deg)
        v_f = transform(v_f, bodies_orientations[target_body.name]['W'], 'z', u.deg)
    return r_f, v_f