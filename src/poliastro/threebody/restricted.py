"""Circular Restricted 3-Body Problem (CR3BP)

    Includes the computation of Lagrange points
"""


import numpy as np
from astropy import units as u

from poliastro._math.optimize import brentq
from poliastro.util import norm


@u.quantity_input(r12=u.km, m1=u.kg, m2=u.kg)
def lagrange_points(r12, m1, m2):
    """Computes the Lagrangian points of CR3BP.

    Computes the Lagrangian points of CR3BP given the distance between two
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
        aux = (1 - pi2) * (xi + pi2) / abs(xi + pi2) ** 3
        aux += pi2 * (xi + pi2 - 1) / abs(xi + pi2 - 1) ** 3
        aux -= xi
        return aux

    lp = np.zeros((5,))

    # L1
    tol = 1e-11  # `brentq` uses a xtol of 2e-12, so it should be covered
    a = -pi2 + tol
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


@u.quantity_input(m1=u.kg, r1=u.km, m2=u.kg, r2=u.km, n=u.one)
def lagrange_points_vec(m1, r1, m2, r2, n):
    """Computes the five Lagrange points in the CR3BP.

    Returns the positions in the same frame of reference as `r1` and `r2`
    for the five Lagrangian points.

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
    assert (
        m1 > m2
    ), "Body 1 is not the main body: it has less mass than the 'secondary' body"

    # Define local frame of reference:
    # Center: main body, NOT the barycenter
    # X-axis: points to the secondary body
    ux = r2 - r1
    r12 = norm(ux)
    ux = ux / r12

    # Y-axis: contained in the orbital plane, perpendicular to x-axis

    def cross(x, y):
        return np.cross(x, y) * x.unit * y.unit

    uy = cross(n, ux)
    uy = uy / norm(uy)

    # Position in x-axis
    x1, x2, x3, x4, x5 = lagrange_points(r12, m1, m2)

    # Position in y-axis
    # L1, L2, L3 are located in the x-axis, so y123 = 0

    # L4 and L5 are points in the plane of rotation which form an equilateral
    # triangle with the two masses (Battin)
    # sqrt(3)/2 = sin(60 deg)
    y4 = np.sqrt(3) / 2 * r12
    y5 = -y4

    # Convert L points coordinates (x,y) to original vectorial base [r1 r2]
    L1 = r1 + ux * x1
    L2 = r1 + ux * x2
    L3 = r1 + ux * x3
    L4 = r1 + ux * x4 + uy * y4
    L5 = r1 + ux * x5 + uy * y5

    return [L1, L2, L3, L4, L5]
