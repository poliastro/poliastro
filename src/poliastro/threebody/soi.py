"""Patched Conics computations

Contains methods to compute interplanetary trajectories approximating the three
body problem with Patched Conics.

"""
from astropy import units as u

from poliastro.constants import J2000
from poliastro.twobody import Orbit


@u.quantity_input(a=u.m)
def laplace_radius(body, a=None):
    """Approximated radius of the Laplace Sphere of Influence (SOI) for a body.

    Parameters
    ----------
    body : `~poliastro.bodies.Body`
           Astronomical body which the SOI's radius is computed for.
    a : float, optional
        Semimajor axis of the body's orbit, default to None (will be computed from ephemerides).

    Returns
    -------
    astropy.units.quantity.Quantity
        Approximated radius of the Sphere of Influence (SOI) [m]

    """
    # Compute semimajor axis at epoch J2000 for the body if it was not
    # introduced by the user
    if a is None:
        try:
            a = Orbit.from_body_ephem(body, J2000).a

        except KeyError:
            raise RuntimeError(
                """To compute the semimajor axis for Moon and Pluto use the JPL ephemeris:

>>> from astropy.coordinates import solar_system_ephemeris
>>> solar_system_ephemeris.set("jpl")"""
            )

    r_SOI = a * (body.k / body.parent.k) ** (2 / 5)

    return r_SOI.decompose()


"""
The laplace radius is focused on finding the body that should be used as the origin when we seek to model the behaviour of the third body in the system.
Laplace radius is computed by the formula:
a\times(\frac{mass_{body}}{mass_{parent}})^{2/5}, where a is the semi-major axis

The hill radius is focused on finding the radius around a body within which a third body will remain in a stable orbit
Hill radius is computed by the formula:
a\times(1-e)\times \sqrt[3]{(\frac{mass_{body}}{3\times mass_{parent}})}, where a is the semi-major axis and e is the eccentricity
"""


@u.quantity_input(a=u.m)
def hill_radius(body, a=None, e=0):
    """Calculates the approximate radius of the Hill sphere of a body.

    Parameters
    ----------
    body : `~poliastro.bodies.Body`
           Astronomical body which the SOI's radius is computed for.
    a : float, optional
        Semimajor axis of the body's orbit, default to None (will be computed from ephemerides).
    e : float, optional
        Eccentricity of the body's orbit, default to 0

    Returns
    -------
    astropy.units.quantity.Quantity
        Approximated radius of the Sphere of Influence (SOI) [m]

    """
    # Compute semimajor axis at epoch J2000 for the body if it was not
    # introduced by the user
    if a is None:
        try:
            a = Orbit.from_body_ephem(body, J2000).a
            # e = Orbit.from_body_ephem(body, J2000).ecc
        except KeyError:
            raise RuntimeError(
                """To compute the semimajor axis and eccentricty for Moon and Pluto use the JPL ephemeris:

>>> from astropy.coordinates import solar_system_ephemeris
>>> solar_system_ephemeris.set("jpl")"""
            )

    r_HR = a * (1 - e) * (body.k / (3 * body.parent.k)) ** (1 / 3)

    return r_HR.decompose()
