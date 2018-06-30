# coding: utf-8

import numpy as np
from scipy.optimize import brentq


def lagrange_points(m1, r1_, m2, r2_, n_, lagrange_point=-1):
    """Function that calculates the five Lagrange points in the
    Restricted Circular 3 Body Problem. Returns the radiovectors in the same
    vectorial base as r1 and r2 for the five Lagrangian points:
    L1, L2, L3, L4, L5

    Parameters
    ----------
    m1 : float
        Mass of the main body. This body is the one with the biggest mass
    r1_ : array
        Radiovector of the main body
    m2 : float
        Mass of the secondary body
    r2_: array
        Radiovector of the secondary body
    n_ : array
        Normal vector to the plane in which the two orbits of the main
        and the secondary body are contained
    lagrange_point: list, optional
        List with the desired Lagrange Points to return. Default value
        is -1, which means all Lagrange Points.
        
    Returns
    -------
    LP:     list of all Lagrange points, unless specified in argument 'lagrange_point'
    """

    m = m2/(m1 + m2)

    r1 = np.asarray(r1_).reshape((3,))
    r2 = np.asarray(r2_).reshape((3,))
    n = np.asarray(n_).reshape((3,))
    # Note: cross product needs the vectors to have this shape

    # Define local reference system:
    # Center: main body
    # x axis: points to the secondary body
    ux = r2 - r1
    r12 = np.linalg.norm(ux)

    # y axis: contained in the orbital plane, perpendicular to x axis
    uy = np.cross(n, ux)

    # Unitary vectors
    ux = ux/r12
    uy = uy/np.linalg.norm(uy)

    print(m)

    if (0.5-m) > 1e-7:

        # Colinear points (y,z always zero in that reference system)
        def eqL1L2L3(x, m): return x - m - (1-m)*x * \
            (abs(x-1)**3) + m*(x-1)*(abs(x)**3)
        # These are obtained by looking for the real roots of eqL1L2L3.

        # L1 is situated between the two main bodies
        x1 = brentq(eqL1L2L3,    0.,   1., args=(m))

        # L2 is situated behind the secondary body (m2,r2)
        x2 = brentq(eqL1L2L3,    1., 1e+7, args=(m))

        # L2 is situated behind the main body (m1,r1)
        x3 = brentq(eqL1L2L3, -1e+7,  -0., args=(m))

        # Triangular points
        # these suppose a problem, since I need to know the plane in which the
        # two main bodies (m1,m2) are orbiting --> n_ clarifies that, as is the
        # normal vector to that plane
        x45 = 0.5
        y45 = np.sqrt(3.)/2.
    else:
        raise ValueError(
            "m = %.5f must be < 0.5, m1 and m2 are too similar or are interchanged" % m)

    # Convert L to original vectors r1 r2 base
    L1 = r1 + ux*r12*x1
    L2 = r1 + ux*r12*x2
    L3 = r1 + ux*r12*x3
    L4 = r1 + ux*r12*x45 + uy*r12*y45
    L5 = r1 + ux*r12*x45 - uy*r12*y45

    LP = [L1, L2, L3, L4, L5]

    if lagrange_point in range(1, 6):
        return LP[lagrange_point-1]
    elif lagrange_point == -1:
        return LP
    else:
        raise ValueError(
            "The value of the argument 'lagrange_point' is not valid.")
