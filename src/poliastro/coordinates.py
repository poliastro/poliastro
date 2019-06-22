"""Functions related to coordinate systems and transformations.

This module complements :py:mod:`astropy.coordinates`.

"""
from math import cos, sin, sqrt

import astropy.units as u
import numpy as np

from poliastro.core.elements import rv2coe
from poliastro.util import norm


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
    DU = u.def_unit("DU", norm(r))
    TU = u.def_unit("TU", np.sqrt((1 * DU) ** 3 / source_body.k))

    r = r.to(DU).value
    v = v.to(DU / TU).value

    p, ecc, inc, _, _, nu = rv2coe(r, v)

    r_pqw = (
        (np.array([cos(nu), sin(nu), 0 * nu]) * p / (1 + ecc * cos(nu))).T * DU
    ).to(u.km)
    v_pqw = ((np.array([-sin(nu), (ecc + cos(nu)), 0]) * sqrt(1 / p)).T * DU / TU).to(
        u.km / u.s
    )

    return r_pqw, v_pqw
