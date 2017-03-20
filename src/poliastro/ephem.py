# coding: utf-8
"""Planetary ephemerides.

"""
from astropy import units as u
from astropy.coordinates import get_body_barycentric_posvel


def get_body_ephem(body, epoch):
    """Position and velocity vectors of a given body at a certain time.

    The vectors are computed with respect to the Solar System barycenter.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    body : str
        Name of the body.
    epoch : astropy.time.Time
        Computation time. Can be scalar or vector.

    Returns
    -------
    r, v : Quantity
        Position and velocity vectors.

    """
    r, v = get_body_barycentric_posvel(body, epoch)
    return r.xyz, v.xyz.to(u.km / u.day)
