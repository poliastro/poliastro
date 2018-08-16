# coding: utf-8

import numpy as np
from scipy.optimize import root


def lagrange_points(r12, m1, m2):
    """Computes the Lagrangian points...

    Parameters
    ----------
    r12 : float
        Collinear distance
    m1 : float
        Mass of the main body
    m2 : float
        Mass of the secondary body

    Returns
    -------
    array(float)
        Distance of the Lagrangian points to the center,
        projected on the collinear line
    """

    # TODO: recheck quintic equations

    rho = r12

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

    # TODO: check boundaries of brentq

    # L1
    rho2 = root(eq_L1, 0)
    rho1 = rho - rho2.x
    l[0] = rho1

    # L2
    rho2 = root(eq_L2, 0)
    rho1 = rho + rho2.x
    l[1] = rho1

    # L3
    # TODO: check -1000 value
    rho1 = root(eq_L3, -1000)
    # rho2 = rho + rho1
    l[2] = - rho1.x

    l[3] = l[4] = 0.5 * r12

    print(l)

    return l


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
    # import matplotlib.pyplot as plt

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

    # plt.figure()
    # plt.scatter(r1[0], r1[1], s=25, marker="o", label="Earth")
    # plt.scatter(r2[0], r2[1], s=25, marker="o", label="Moon")
    # plt.scatter(x, y, marker="x", label="L points")
    # plt.legend(loc="best")
    # plt.show()
