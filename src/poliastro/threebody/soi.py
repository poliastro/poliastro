"""Sphere of Influence

Contains methods to compute radius of sphere of influence.

Laplace Radius: A laplace sphere of influence (SOI) in astrodynamics and
astronomy is the oblate-spheroid-shaped region around a celestial body where
the primary gravitational influence on an orbiting object is that body.  This
is usually used to describe the areas in the Solar System where planets
dominate the orbits of surrounding objects such as moons, despite the presence
of the much more massive but distant Sun. In the patched conic approximation,
used in estimating the trajectories of bodies moving between the neighbourhoods
of different masses using a two body approximation, ellipses and hyperbolae,
the laplace radius is taken as the boundary where the trajectory switches which
mass field it is influenced by.  The result is:

.. math::
    a\\left(\\frac{m}{M}\\right)^{\\frac{2}{5}}

Hill Radius: In the three body problem, if that third object stays within an
extremely complex boundary called the Roche lobe, the orbit of that third
object about the smaller body will be stable for at least some amount of time.
The Roche lobe just touches the L1 and L2 points and fans out from there.
George Hill used the L1 point to define a sphere that approximated the Roche
lobe. This is still intractable; the L1 point is defined by a fifth order
polynomial that cannot be solved in terms of the elementary functions. Hill
further simplified things by realizing that a simple cubic equation yields a
very good approximation of that intractable fifth order equation.  The result
is:

.. math::
    a\\left(\\frac{m}{3M}\\right)^{\\frac{1}{3}}

"""
from poliastro.twobody.mean_elements import get_mean_elements


def laplace_radius(body):
    """Approximated radius of the Laplace Sphere of Influence (SOI) for a body.

    Parameters
    ----------
    body : `~poliastro.bodies.Body`
        Astronomical body which the SOI's radius is computed for.

    Returns
    -------
    astropy.units.quantity.Quantity
        Approximated radius of the Laplace Sphere of Influence

    """
    a = get_mean_elements(body).a
    r_SOI = a * (body.k / body.parent.k) ** (2 / 5)

    return r_SOI.decompose()


def hill_radius(body):
    """Approximated radius of the Hill Sphere of Influence (SOI) for a body.

    Parameters
    ----------
    body : `~poliastro.bodies.Body`
        Astronomical body which the SOI's radius is computed for.

    Returns
    -------
    astropy.units.quantity.Quantity
        Approximated radius of the Hill Sphere of Influence

    """
    mean_elements = get_mean_elements(body)
    a = mean_elements.a
    ecc = mean_elements.ecc

    mass_ratio = body.k / (3 * body.parent.k)
    r_SOI = a * (1 - ecc) * (mass_ratio ** (1 / 3))

    return r_SOI.decompose()
