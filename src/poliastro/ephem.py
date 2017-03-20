# coding: utf-8
"""Planetary ephemerides.

"""
from astropy import units as u
from astropy.coordinates import get_body_barycentric_posvel


# Stars
SUN = "sun"

# Planets
MERCURY = "mercury"
VENUS = "venus"
EARTH = "earth"
MARS = "mars"
JUPITER = "jupiter"
SATURN = "saturn"
URANUS = "uranus"
NEPTUNE = "neptune"

# Other bodies
PLUTO = "pluto"


def get_body_ephem(body, epoch):
    """Position and velocity vectors of a given body at a certain time.

    The vectors are computed with respect to the Solar System barycenter.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    body : Body, str or int
        Planetary body. Can be a Body object, an string or an integer.
    epoch : astropy.time.Time
        Computation time. Can be scalar or vector.

    Returns
    -------
    r, v : Quantity
        Position and velocity vectors.

    """
    try:
        body_str = body.name
    except AttributeError:  # The body is not a `Body` object, assume str
        body_str = body

    r, v = get_body_barycentric_posvel(body_str, epoch)
    return r.xyz, v.xyz.to(u.km / u.day)
