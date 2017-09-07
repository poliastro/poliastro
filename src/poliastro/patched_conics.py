# coding: utf-8
"""Patched Conics Computations

Contains methods to compute interplanetary trajectories approximating the three
body problem with Patched Conics.
"""

from astropy import units as u
from poliastro.twobody import Orbit
from poliastro.constants import J2000


@u.quantity_input(a=u.m)
def compute_soi(body, a=None):
    """Approximated radius of the Laplace Sphere of Influence (SOI) for a body.

    Parameters
    ----------
    body : `~poliastro.bodies.Body`
           Astronomical body which the SOI's radius is computed for

    a : float or None, optional
        Semimajor Axis of the body's orbit

    Returns
    -------
    astropy.units.quantity.Quantity
        Approximated radius of the Sphere of Influence (SOI) [m]
    """
    # Compute semimajor axis at epoch J2000 for the body if it was not
    # introduced by the user
    if a is None:
        ss = Orbit.from_body_ephem(body, J2000)
        a = ss.a

    r_SOI = a * (body.k / body.parent.k)**(2 / 5)

    return r_SOI.decompose()
