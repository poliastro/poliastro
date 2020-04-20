"""Functions related to coordinate systems and transformations.

This module complements :py:mod:`astropy.coordinates`.

"""
from math import cos, sin, sqrt

import astropy.units as u
import numpy as np

from poliastro.core.elements import rv2coe


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
