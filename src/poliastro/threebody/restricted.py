"""Restricted Circular 3-Body Problem (RC3BP)

    Includes the computation of the Lagrange points
"""


import numpy as np

from scipy.optimize import root, brentq

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

    l = np.zeros((5,))

    # L1
    xi = root(eq_L123, 1 - 1.9*pi2).x
    l[0] = xi + pi2

    # L2
    xi = brentq(eq_L123, 1, 1.5)
    l[1] = xi + pi2

    # L3
    xi = brentq(eq_L123, -1.5, -1)
    l[2] = xi + pi2

    # L4, L5
    # L4 and L5 are points in the plane of rotation which form an equilateral
    # triangle with the two masses (Battin)
    # (0.5 = cos(60 deg))
    l[3] = l[4] = 0.5

    return l * r12


def lagrange_points_battin(r12, m1, m2):
    """Computes the Lagrangian points of RC3BP given the distance between two
    bodies and their masses.

    It uses the formulation found in Problem 8-5 of Battin, Richard H. 'An 
    introduction to the mathematics and methods of astrodynamics', AIAA, 1999.

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

    rho = r12.value
    m1 = m1.value
    m2 = m2.value

    def eq_L1(rho2):
        aux = 3 * rho**2 * rho2 - 3 * rho * rho2**2 + rho2**3
        aux *= rho2**2
        aux /= rho**3 - rho2**3
        aux /= (rho - rho2)**2
        aux -= m2 / m1
        return aux

    def eq_L2(rho2):
        aux = 3 * rho**2 * rho2 + 3 * rho * rho2**2 + rho2**3
        aux *= rho2**2
        aux /= rho**3 - rho2**3
        aux /= (rho + rho2)**2
        aux -= m2 / m1
        return aux

    def eq_L3(rho1):
        aux = 3 * rho**2 * rho1 + 3 * rho * rho1**2 + rho1**3
        aux *= rho1**2
        aux /= rho**3 - rho1**3
        aux /= (rho + rho1)**2
        aux -= m1 / m2
        return aux

    l = np.zeros((5,))

    # L1
    rho2 = root(eq_L1, 0).x
    rho1 = rho - rho2
    l[0] = rho1

    # L2
    rho2 = root(eq_L2, 0).x
    rho1 = rho + rho2
    l[1] = rho1

    # L3
    rho1 = root(eq_L3, -1000).x
    # rho2 = rho + rho1
    l[2] = - rho1

    # TODO: add checks to the results returned by `root`
    # TODO: reassure that the initial values given to `root` are good for all cases

    # L4, L5
    # L4 and L5 are points in the plane of rotation which form an equilateral
    # triangle with the two masses (Battin)
    # (0.5 = cos(60 deg))
    l[3] = l[4] = 0.5 * rho

    return l * r12.unit


def lagrange_points_vec(m1, r1, m2, r2, n):
    """Computes the five Lagrange points in the RC3BP. Returns the positions
    in the same vectorial base as `r1` and `r2` for the five Lagrangian points.

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
    assert m1 > m2, "Body 1 is not the main body: it has less mass that the 'secondary' body"

    # Define local reference system:
    # Center: main body, NOT the barycenter
    # x axis: points to the secondary body
    ux = r2 - r1
    r12 = norm(ux)
    ux = ux / r12

    # y axis: contained in the orbital plane, perpendicular to x axis

    def cross(x, y):
        return np.cross(x, y) * x.unit * y.unit

    uy = cross(n, ux)
    uy = uy / norm(uy)

    # position in x-axis
    x1, x2, x3, x4, x5 = lagrange_points(r12, m1, m2)

    # position in y axis
    # L1, L2, L3 are located in the x-axis, so y123 = 0

    # L4 and L5 are points in the plane of rotation which form an equilateral
    # triangle with the two masses (Battin)
    # sqrt(3)/2 = sin(60 deg)
    y4 = np.sqrt(3) / 2 * r12
    y5 = - y4

    # Convert L points (x,y) to original vectors [r1 r2] base
    L1 = r1 + ux * x1
    L2 = r1 + ux * x2
    L3 = r1 + ux * x3
    L4 = r1 + ux * x4 + uy * y4
    L5 = r1 + ux * x5 + uy * y5

    return [L1, L2, L3, L4, L5]


if __name__ == "__main__":

    from astropy import units as u
    from astropy.constants import G
    from poliastro.constants import GM_earth, GM_moon

    # ORIGIN = "barycenter"
    ORIGIN = "main body"

    # Distance Earth - Moon
    r12 = 384400

    # Earth
    m1 = GM_earth / G

    # Moon
    m2 = GM_moon / G

    if ORIGIN == "barycenter":
        x1 = - r12 * m2 / (m1 + m2)
        x2 = r12 + x1
    elif ORIGIN == "main body":
        x1 = 0
        x2 = r12

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
