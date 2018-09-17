"""Restricted Circular 3-Body Problem (RC3BP)

    Includes the computation of:
    * Lagrange points
"""


import numpy as np

from scipy.optimize import brentq

from poliastro.util import norm


def lagrange_points(r12, m1, m2):
    """Computes the Lagrangian points of RC3BP given the distance between two
    bodies and their masses.

    It uses the formulation found in Eq. (2.204) of Curtis, Howard. 'Orbital
    mechanics for engineering students'. Elsevier, 3rd Edition.

    Parameters
    ----------
    r12 : ~astropy.units.Quantity
        Distance between the two bodies
    m1 : ~astropy.units.Quantity
        Mass of the main body
    m2 : ~astropy.units.Quantity
        Mass of the secondary body

    Returns
    -------
    ~astropy.units.Quantity
        Distance of the Lagrangian points to the main body,
        projected on the axis main body - secondary body
    """

    pi2 = (m2 / (m1 + m2)).value

    def eq_L123(xi):
        aux = (1 - pi2) * (xi + pi2) / abs(xi + pi2)**3
        aux += pi2 * (xi + pi2 - 1) / abs(xi + pi2 - 1)**3
        aux -= xi
        return aux

    lp = np.zeros((5,))

    # L1
    tol = 1e-11  # `brentq` uses a xtol of 2e-12, so it should be covered
    a = - pi2 + tol
    b = 1 - pi2 - tol
    xi = brentq(eq_L123, a, b)
    lp[0] = xi + pi2

    # L2
    xi = brentq(eq_L123, 1, 1.5)
    lp[1] = xi + pi2

    # L3
    xi = brentq(eq_L123, -1.5, -1)
    lp[2] = xi + pi2

    # L4, L5
    # L4 and L5 are points in the plane of rotation which form an equilateral
    # triangle with the two masses (Battin)
    # (0.5 = cos(60 deg))
    lp[3] = lp[4] = 0.5

    return lp * r12


def lagrange_points_vec(m1, r1, m2, r2, n):
    """Computes the five Lagrange points in the RC3BP. Returns the positions
    in the same frame of reference as `r1` and `r2` for the five Lagrangian
    points.

    Parameters
    ----------
    m1 : ~astropy.units.Quantity
        Mass of the main body. This body is the one with the biggest mass.
    r1 : ~astropy.units.Quantity
        Position of the main body.
    m2 : ~astropy.units.Quantity
        Mass of the secondary body.
    r2 : ~astropy.units.Quantity
        Position of the secondary body.
    n : ~astropy.units.Quantity
        Normal vector to the orbital plane.

    Returns
    -------
    list:
        Position of the Lagrange points: [L1, L2, L3, L4, L5]
        The positions are of type ~astropy.units.Quantity
    """

    # Check Body 1 is the main body
    assert m1 > m2, "Body 1 is not the main body: it has less mass than the 'secondary' body"

    # Define local frame of reference:
    # Center: main body, NOT the barycenter
    # x-axis: points to the secondary body
    ux = r2 - r1
    r12 = norm(ux)
    ux = ux / r12

    # y-axis: contained in the orbital plane, perpendicular to x-axis

    def cross(x, y):
        return np.cross(x, y) * x.unit * y.unit

    uy = cross(n, ux)
    uy = uy / norm(uy)

    # position in x-axis
    x1, x2, x3, x4, x5 = lagrange_points(r12, m1, m2)

    # position in y-axis
    # L1, L2, L3 are located in the x-axis, so y123 = 0

    # L4 and L5 are points in the plane of rotation which form an equilateral
    # triangle with the two masses (Battin)
    # sqrt(3)/2 = sin(60 deg)
    y4 = np.sqrt(3) / 2 * r12
    y5 = - y4

    # Convert L points coordinates (x,y) to original vectorial base [r1 r2]
    L1 = r1 + ux * x1
    L2 = r1 + ux * x2
    L3 = r1 + ux * x3
    L4 = r1 + ux * x4 + uy * y4
    L5 = r1 + ux * x5 + uy * y5

    return [L1, L2, L3, L4, L5]


def lagrange_points_from_body(body):
    """Computes the five Lagrange points given a `Body`.

    Parameters
    ----------
    body : Body
        Orbiting body. 
        Example: For Earth-Moon, pass Moon, because the attractor body is already known.

    Returns
    -------
    list
        Position of the Lagrange points: [L1, L2, L3, L4, L5]
        The positions are of type ~astropy.units.Quantity
    """
    from poliastro.twobody import Orbit
    from astropy import units as u

    body_orbit = Orbit.from_body_ephem(body)
    body_orbit_parent_r = np.zeros(3) * u.km  # origin in main body

    return lagrange_points_vec(m1=body.parent.mass,
                               r1=body_orbit_parent_r,
                               m2=body.mass,
                               r2=body_orbit.r,
                               n=body_orbit.h_vec)


def test_sun_earth():
    from astropy import units as u
    from astropy.constants import G, au
    from poliastro.constants import GM_earth, GM_sun

    # ORIGIN = "barycenter"
    ORIGIN = "main body"

    # Distance Sun - Earth
    r12 = au.to(u.km).value

    # Sun (1)
    m1 = GM_sun / G

    # Earth (2)
    m2 = GM_earth / G

    if ORIGIN == "barycenter":
        x1 = - r12 * m2 / (m1 + m2)
        x2 = r12 + x1
    elif ORIGIN == "main body":
        x1 = 0.
        x2 = r12

    # Positions
    r1 = np.array([x1, 0, 0]) * u.km
    r2 = np.array([x2, 0, 0]) * u.km

    # normal vector
    n = np.array([0., 0, 1]) * u.one

    lp = lagrange_points_vec(m1, r1, m2, r2, n)

    for p in lp:
        print("{:+8.0f} {:+8.0f} {:+8.0f}".format(p[0], p[1], p[2]))

    x = [p[0] for p in lp]
    y = [p[1] for p in lp]

    # Figure
    import matplotlib.pyplot as plt
    plt.figure()
    plt.scatter(0, 0, marker="+", c="k", label="Origin")
    plt.scatter(r1[0], r1[1], s=100, marker="$" + u"\u2609" + "$",
                label="Sun", c="k")
    plt.scatter(r2[0], r2[1], s=100, marker="$" + u"\u2641" + "$",
                label="Earth", c="k")
    for i in range(0, 5):
        plt.scatter(x[i], y[i], marker="$L" + "{:d}$".format(i + 1),
                    c="k", s=100)
    plt.legend(loc="best")
    plt.title("Sun-Earth Lagrangian points")
    plt.xlabel("[km]")
    plt.ylabel("[km]")
    plt.tight_layout()
    plt.show()
    plt.close('all')


def test_earth_moon():
    from astropy import units as u
    from astropy.constants import G
    from poliastro.constants import GM_earth, GM_moon

    # ORIGIN = "barycenter"
    ORIGIN = "main body"

    # Distance Earth - Moon
    r12 = 384400

    # Earth (1)
    m1 = GM_earth / G

    # Moon (2)
    m2 = GM_moon / G

    if ORIGIN == "barycenter":
        x1 = - r12 * m2 / (m1 + m2)
        x2 = r12 + x1
    elif ORIGIN == "main body":
        x1 = 0
        x2 = r12

    # Positions
    r1 = np.array([x1, 0, 0]) * u.km
    r2 = np.array([x2, 0, 0]) * u.km

    # normal vector
    n = np.array([0., 0, 1]) * u.one

    lp = lagrange_points_vec(m1, r1, m2, r2, n)

    for p in lp:
        print("{:+8.0f} {:+8.0f} {:+8.0f}".format(p[0], p[1], p[2]))

    x = [p[0] for p in lp]
    y = [p[1] for p in lp]

    # Figure
    import matplotlib.pyplot as plt
    plt.figure()
    plt.scatter(0, 0, marker="+", c="k", label="Origin")
    plt.scatter(r1[0], r1[1], s=100, marker="$" + u"\u2641" + "$",
                label="Earth", c="k")
    plt.scatter(r2[0], r2[1], s=100, marker="$" + u"\u263E" + "$",
                label="Moon", c="k")
    for i in range(0, 5):
        plt.scatter(x[i], y[i], marker="$L" + "{:d}$".format(i + 1),
                    c="k", s=100)
    plt.legend(loc="best")
    plt.title("Earth-Moon Lagrangian points")
    plt.xlabel("[km]")
    plt.ylabel("[km]")
    plt.tight_layout()
    plt.show()
    plt.close('all')
