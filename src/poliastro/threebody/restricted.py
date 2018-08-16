# coding: utf-8

import numpy as np
from scipy.optimize import brentq


def lagrange_points(r12, m1, m2):
    """Computes the Lagrangian points proe

    Parameters
    ----------
    r12 : float
        Collinear distance
    m1 : float
        Mass of the main body
    m2 : float
        Mass of the secondary body

    Raises
    ------
    ValueError
        If the ratio $m_2 / (m_1 + m_2)$ is less than 0.5

    Returns
    -------
    array(float)
        Distance of the Lagrangian points to the center,
        projected on the collinear line
    """

    m = m2 / (m1 + m2)

    if m < 0.5:

        def eq_collinear(x):
            r = x - m - (1 - m) * x * (abs(x - 1) ** 3)
            r += m * (x - 1) * (abs(x)**3)
            return r

        l = np.zeros((5,))

        # L1 is situated between the two main bodies
        l[0] = brentq(eq_collinear, 0., 1.)

        # L2 is situated behind the secondary body (m2,r2)
        l[1] = brentq(eq_collinear, 1., 1e+7)

        # L3 is situated behind the main body (m1,r1)
        l[2] = brentq(eq_collinear, -1e+7, -0.)

        l[3] = l[4] = 0.5

        return l * r12

    else:
        raise ValueError(
            "m = {:.5f} must be < 0.5, m1 and m2 are ".format(m) +
            "too similar or are interchanged")


def lagrange_points_vec(m1, r1_, m2, r2_, n_):
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

    Returns
    -------
    LP: list of all Lagrange points
    """

    r1 = np.asarray(r1_).reshape((3,))
    r2 = np.asarray(r2_).reshape((3,))
    n = np.asarray(n_).reshape((3,))

    # Define local reference system:
    # Center: main body
    # x axis: points to the secondary body
    ux = r2 - r1
    r12 = np.linalg.norm(ux)

    # y axis: contained in the orbital plane, perpendicular to x axis
    uy = np.cross(n, ux)

    # Unitary vectors
    ux = ux / r12
    uy = uy / np.linalg.norm(uy)

    x1, x2, x3, x4, x5 = lagrange_points(r12, m1, m2)

    y45 = np.sqrt(3) / 2

    # Convert L to original vectors r1 r2 base
    L1 = r1 + ux * x1
    L2 = r1 + ux * x2
    L3 = r1 + ux * x3
    L4 = r1 + ux * x4 + uy * y45
    L5 = r1 + ux * x5 - uy * y45

    return [L1, L2, L3, L4, L5]


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    # Earth
    r1 = np.array([0, 0, 0])
    m1 = 5.972e24

    # Moon
    m2 = 7.348e22
    d = 384400
    r2 = np.array([d, 0, 0])

    n = np.array([0., 0, 1])

    res = lagrange_points_vec(m1, r1, m2, r2, n)
    res1 = np.array(res) / d

    [print("{:+10.5e} , {:+10.5e} , {:+10.5e}".format(r[0], r[1], r[2]))
     for r in res1]

    x = np.array([r[0] for r in res])
    y = np.array([r[1] for r in res])

    plt.figure()
    plt.scatter(r1[0], r1[1], s=25, marker="o", label="Earth")
    plt.scatter(r2[0], r2[1], s=25, marker="o", label="Moon")
    plt.scatter(x, y, marker="x", label="L points")
    plt.legend(loc="best")
    plt.show()
