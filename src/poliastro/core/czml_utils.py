import numpy as np
from numba import njit as jit

from poliastro._math.linalg import norm


@jit
def intersection_ellipsoid_line(x, y, z, u1, u2, u3, a, b, c):
    """Intersection of an ellipsoid defined by its axes a, b, c with the
    line p + λu.

    Parameters
    ----------
    x, y, z : float
        A point of the line
    u1, u2, u3 : float
        The line vector
    a, b, c : float
        The ellipsoidal axises

    Returns
    -------
    p0, p1: numpy.ndarray
        This returns both of the points intersecting the ellipsoid.

    """
    # Get rid of one parameter by translating the line's direction vector
    k, m = u2 / u1, u3 / u1
    t0 = (
        -(a**2) * b**2 * m * z
        - a**2 * c**2 * k * y
        - b**2 * c**2 * x
        + np.sqrt(
            a**2
            * b**2
            * c**2
            * (
                a**2 * b**2 * m**2
                + a**2 * c**2 * k**2
                - a**2 * k**2 * z**2
                + 2 * a**2 * k * m * y * z
                - a**2 * m**2 * y**2
                + b**2 * c**2
                - b**2 * m**2 * x**2
                + 2 * b**2 * m * x * z
                - b**2 * z**2
                - c**2 * k**2 * x**2
                + 2 * c**2 * k * x * y
                - c**2 * y**2
            )
        )
    ) / (a**2 * b**2 * m**2 + a**2 * c**2 * k**2 + b**2 * c**2)
    t1 = (
        a**2 * b**2 * m * z
        + a**2 * c**2 * k * y
        + b**2 * c**2 * x
        + np.sqrt(
            a**2
            * b**2
            * c**2
            * (
                a**2 * b**2 * m**2
                + a**2 * c**2 * k**2
                - a**2 * k**2 * z**2
                + 2 * a**2 * k * m * y * z
                - a**2 * m**2 * y**2
                + b**2 * c**2
                - b**2 * m**2 * x**2
                + 2 * b**2 * m * x * z
                - b**2 * z**2
                - c**2 * k**2 * x**2
                + 2 * c**2 * k * x * y
                - c**2 * y**2
            )
        )
    ) / (a**2 * b**2 * m**2 + a**2 * c**2 * k**2 + b**2 * c**2)
    p0, p1 = np.array([x + t0, y + k * t0, z + m * t0]), np.array(
        [x - t1, y - t1 * k, z - t1 * m]
    )

    return p0, p1


@jit
def project_point_on_ellipsoid(x, y, z, a, b, c):
    """Return the projection of a point on an ellipsoid.

    Parameters
    ----------
    x, y, z : float
        Cartesian coordinates of point
    a, b, c : float
        Semi-axes of the ellipsoid

    """
    p1, p2 = intersection_ellipsoid_line(x, y, z, x, y, z, a, b, c)

    norm_1 = norm(np.array([p1[0] - x, p1[1] - y, p1[2] - z]))
    norm_2 = norm(np.array([p2[0] - x, p2[1] - y, p2[2] - z]))

    return p1 if norm_1 <= norm_2 else p2
