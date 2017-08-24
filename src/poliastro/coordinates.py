"""Functions related to coordinate systems and transformations.

This module complements :py:mod:`astropy.coordinates`.

"""

import math

import astropy.time as time
import astropy.coordinates as coordinates
import astropy.units as u

from poliastro.util import transform
from poliastro.constants import J2000


def mercury_rot_elements(epoch=J2000):
    """Provides Mercury rotational elements given an epoch.

    Provides north pole of body and angle to prime meridian

    Parameters
    ----------
    epoch : ~astropy.time.Time, optional
        Epoch, default to J2000.

    Returns
    -------
    ra, dec, W: tuple (~astropy.units.Quantity)
        Right ascension and declination of north pole, and angle of the prime meridian.

    """

    T = (epoch.tdb - J2000).to('day').value / 36525
    d = (epoch.tdb - J2000).to('day').value

    ra = (281.01 - 0.033 * T) * u.deg
    dec = (61.45 - 0.005 * T) * u.deg
    W = (329.548 + 6.1385025 * d) * u.deg

    return ra, dec, W


def venus_rot_elements(epoch=J2000):
    """Provides Venus rotational elements given an epoch.

    Provides north pole of body and angle to prime meridian

    Parameters
    ----------
    epoch : ~astropy.time.Time, optional
        Epoch, default to J2000.

    Returns
    -------
    ra, dec, W: tuple (~astropy.units.Quantity)
        Right ascension and declination of north pole, and angle of the prime meridian.

    """

    d = (epoch.tdb - J2000).to('day').value

    ra = 272.76 * u.deg
    dec = 67.16 * u.deg
    W = (160.20 - 1.4813688 * d) * u.deg

    return ra, dec, W


def earth_rot_elements(epoch=J2000):
    """Provides Earth rotational elements given an epoch.

    Provides north pole of body and angle to prime meridian

    Parameters
    ----------
    epoch : ~astropy.time.Time, optional
        Epoch, default to J2000.

    Returns
    -------
    ra, dec, W: tuple (~astropy.units.Quantity)
        Right ascension and declination of north pole, and angle of the prime meridian.

    """
    T = (epoch.tdb - J2000).to('day').value / 36525
    d = (epoch.tdb - J2000).to('day').value

    ra = (0.00 - 0.641 * T) * u.deg
    dec = (90.00 - 0.557 * T) * u.deg
    W = (190.147 + 360.9856235 * d) * u.deg

    return ra, dec, W


def mars_rot_elements(epoch=J2000):
    """Provides Mars rotational elements given an epoch.

    Provides north pole of body and angle to prime meridian

    Parameters
    ----------
    epoch : ~astropy.time.Time, optional
        Epoch, default to J2000.

    Returns
    -------
    ra, dec, W: tuple (~astropy.units.Quantity)
        Right ascension and declination of north pole, and angle of the prime meridian.

    """
    T = (epoch.tdb - J2000).to('day').value / 36525
    d = (epoch.tdb - J2000).to('day').value

    ra = (317.68143 - 0.1061 * T) * u.deg
    dec = (52.88650 - 0.0609 * T) * u.deg
    W = (176.630 + 350.89198226 * d) * u.deg

    return ra, dec, W


def jupiter_rot_elements(epoch=J2000):
    """Provides Jupiter rotational elements given an epoch.

    Provides north pole of body and angle to prime meridian

    Parameters
    ----------
    epoch : ~astropy.time.Time, optional
        Epoch, default to J2000.

    Returns
    -------
    ra, dec, W: tuple (~astropy.units.Quantity)
        Right ascension and declination of north pole, and angle of the prime meridian.

    """
    T = (epoch.tdb - J2000).to('day').value / 36525
    d = (epoch.tdb - J2000).to('day').value

    Ja = (99.360714 + 4850.4046 * T) * u.deg
    Jb = (175.895369 + 1191.9605 * T) * u.deg
    Jc = (300.323162 + 262.5475 * T) * u.deg
    Jd = (114.012305 + 6070.2476 * T) * u.deg
    Je = (49.511251 + 64.3000 * T) * u.deg

    ra = (268.056595 - 0.006499 * T + 0.000117 * math.sin(Ja) + 0.000938 * math.sin(Jb) +
          0.001432 * math.sin(Jc) + 0.000030 * math.sin(Jd) + 0.002150 * math.sin(Je)) * u.deg
    dec = (64.495303 + 0.002413 * T + 0.000050 * math.cos(Ja) + 0.000404 * math.cos(Jb) +
           0.000617 * math.cos(Jc) - 0.000013 * math.cos(Jd) + 0.000926 * math.cos(Je)) * u.deg
    W = (284.95 + 870.5366420 * d) * u.deg

    return ra, dec, W


def saturn_rot_elements(epoch=J2000):
    """Provides Saturn rotational elements given an epoch.

    Provides north pole of body and angle to prime meridian

    Parameters
    ----------
    epoch : ~astropy.time.Time, optional
        Epoch, default to J2000.

    Returns
    -------
    ra, dec, W: tuple (~astropy.units.Quantity)
        Right ascension and declination of north pole, and angle of the prime meridian.

    """
    T = (epoch.tdb - J2000).to('day').value / 36525
    d = (epoch.tdb - J2000).to('day').value

    ra = (40.589 - 0.036 * T) * u.deg
    dec = (83.537 - 0.004 * T) * u.deg
    W = (38.90 + 810.7939024 * d) * u.deg

    return ra, dec, W


def uranus_rot_elements(epoch=J2000):
    """Provides Uranus rotational elements given an epoch.

    Provides north pole of body and angle to prime meridian

    Parameters
    ----------
    epoch : ~astropy.time.Time, optional
        Epoch, default to J2000.

    Returns
    -------
    ra, dec, W: tuple (~astropy.units.Quantity)
        Right ascension and declination of north pole, and angle of the prime meridian.

    """
    d = (epoch.tdb - J2000).to('day').value

    ra = 257.311 * u.deg
    dec = -15.175 * u.deg
    W = (203.81 - 501.1600928 * d) * u.deg

    return ra, dec, W


def neptune_rot_elements(epoch=J2000):
    """Provides Uranus rotational elements given an epoch.

    Provides north pole of body and angle to prime meridian

    Parameters
    ----------
    epoch : ~astropy.time.Time, optional
        Epoch, default to J2000.

    Returns
    -------
    ra, dec, W: tuple (~astropy.units.Quantity)
        Right ascension and declination of north pole, and angle of the prime meridian.

    """
    T = (epoch.tdb - J2000).to('day').value / 36525
    N = 357.85 + 52.316 * T
    d = (epoch.tdb - J2000).to('day').value

    ra = (299.36 + 0.70 * math.sin(N)) * u.deg
    dec = (43.46 - 0.51 * math.cos(N)) * u.deg
    W = (253.18 + 536.3128492 * d - 0.48 * math.sin(N)) * u.deg

    return ra, dec, W


def pluto_rot_elements(epoch=J2000):
    """Provides Pluto rotational elements given an epoch.

    Provides north pole of body and angle to prime meridian

    Parameters
    ----------
    epoch : ~astropy.time.Time, optional
        Epoch, default to J2000.

    Returns
    -------
    ra, dec, W: tuple (~astropy.units.Quantity)
        Right ascension and declination of north pole, and angle of the prime meridian.

    """
    d = (epoch.tdb - J2000).to('day').value

    ra = 312.993 * u.deg
    dec = 6.163 * u.deg
    W = (190.147 + 360.9856235 * d) * u.deg

    return ra, dec, W


def rot_elements_from_body(body, epoch=J2000):
    """Provides rotational elements of a given body at epoch.

    Provides north pole of body and angle to prime meridian

    Parameters
    ----------
    body : ~poliastro.bodies.Body
        Body.
    epoch : ~astropy.time.Time, optional
        Epoch, default to J2000.

    Returns
    -------
    ra, dec, W: tuple (~astropy.units.Quantity)
        Right ascension and declination of north pole, and angle of the prime meridian.

    """

    options = {'Mercury': mercury_rot_elements,
               'Venus': venus_rot_elements,
               'Earth': earth_rot_elements,
               'Mars': mars_rot_elements,
               'Jupiter': jupiter_rot_elements,
               'Saturn': saturn_rot_elements,
               'Uranus': uranus_rot_elements,
               'Neptune': neptune_rot_elements,
               'Pluto': pluto_rot_elements,
               }

    return options[body.name](epoch)


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

    ra, dec, W = rot_elements_from_body(source_body, epoch)
    if rotate_meridian:
        r = transform(r, -W.value, 'z', W.unit)
        v = transform(v, -W.value, 'z', W.unit)

    r_trans1 = transform(r, -(90 - dec.value), 'x', dec.unit)
    r_trans2 = transform(r_trans1, -(90 + ra.value), 'z', ra.unit)

    v_trans1 = transform(v, -(90 - dec.value), 'x', dec.unit)
    v_trans2 = transform(v_trans1, -(90 + ra.value), 'z', ra.unit)

    icrs_frame_pos_coord, icrs_frame_vel_coord = coordinates.get_body_barycentric_posvel(source_body.name, time=epoch)

    r_f = icrs_frame_pos_coord.xyz + r_trans2
    v_f = icrs_frame_vel_coord.xyz + v_trans2

    # icrs_coord = coordinates.ICRS(x=r_f[0], y=r_f[1], z=r_f[2], v_x=v_f[0], v_y=v_f[1],
    #                              v_z=v_f[2], representation='cartesian')
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

    ra, dec, W = rot_elements_from_body(target_body, epoch)

    icrs_frame_pos_coord, icrs_frame_vel_coord = coordinates.get_body_barycentric_posvel(target_body.name, time=epoch)

    r_trans1 = r - icrs_frame_pos_coord.xyz
    r_trans2 = transform(r_trans1, (90 + ra.value), 'z', ra.unit)
    r_f = transform(r_trans2, (90 - dec.value), 'x', dec.unit)

    v_trans1 = v - icrs_frame_vel_coord.xyz
    v_trans2 = transform(v_trans1, (90 + ra.value), 'z', ra.unit)
    v_f = transform(v_trans2, (90 - dec.value), 'x', dec.unit)

    if rotate_meridian:
        r_f = transform(r_f, W.value, 'z', W.unit)
        v_f = transform(v_f, W.value, 'z', W.unit)

    return r_f.to(r.unit), v_f.to(v.unit)
