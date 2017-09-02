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
    """Computes the approximated radius of the Laplace Sphere of Influence (SOI)
 for a body (:py:class:`~Body` class).
    """
    # Compute semimajor axis at epoch J2000 for the body if it was not
    # introduced by the user
    if a == None:
        ss = Orbit.from_body_ephem(body, J2000)
        a = ss.a

    try:
        r_soi = a * (body.k / body.parent.k)**(2 / 5)
        return r_soi.decompose()
    except AttributeError:
        print('Body has no parent body assigned.')
