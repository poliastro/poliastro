import numpy as np


def ellipsoidal_to_cartesian(a, b, lat, lon):
    """

    Converts ellipsoidal coordinates (assuming zero elevation)
    to the cartesian.

    Parameters
    ----------
    a: float
        semi-major axis
    b: float
        semi-minor axis
    lat: float
        latitude
    lon: float
        longtitude

    Returns
    -------
    3d Cartesian coordinates in the form (x, y, z)
    """

    # radius of curvature
    cr = a / np.sqrt(1 - (1 - (a / b) ** 2) * np.sin(lat) ** 2)

    x = cr * np.cos(lat) * np.cos(lon)
    y = cr * np.cos(lat) * np.sin(lon)
    z = (a / b) ** 2 * cr * np.sin(lat)
    return x, y, z


def intersection_ellipsoid_line(x, y, z, u1, u2, u3, a, b, c):
    """
    Intersection of an ellipsoid defined by its axises a, b, c with the
    line p + λu
    Parameters
    ----------
    x, y, z: float
        a point of the line
    u1, u2, u3: float
        the line vector
    a, b, c: float
        the ellipsoidal axises
    Returns
    -------
    This returns both of the points intersecting the ellipsoid.
    """
    # Get rid of one parameter by translating the line's direction vector
    k, m = u2 / u1, u3 / u1
    t0 = (
        -a ** 2 * b ** 2 * m * z
        - a ** 2 * c ** 2 * k * y
        - b ** 2 * c ** 2 * x
        + np.sqrt(
            a ** 2
            * b ** 2
            * c ** 2
            * (
                a ** 2 * b ** 2 * m ** 2
                + a ** 2 * c ** 2 * k ** 2
                - a ** 2 * k ** 2 * z ** 2
                + 2 * a ** 2 * k * m * y * z
                - a ** 2 * m ** 2 * y ** 2
                + b ** 2 * c ** 2
                - b ** 2 * m ** 2 * x ** 2
                + 2 * b ** 2 * m * x * z
                - b ** 2 * z ** 2
                - c ** 2 * k ** 2 * x ** 2
                + 2 * c ** 2 * k * x * y
                - c ** 2 * y ** 2
            )
        )
    ) / (a ** 2 * b ** 2 * m ** 2 + a ** 2 * c ** 2 * k ** 2 + b ** 2 * c ** 2)
    t1 = (
        a ** 2 * b ** 2 * m * z
        + a ** 2 * c ** 2 * k * y
        + b ** 2 * c ** 2 * x
        + np.sqrt(
            a ** 2
            * b ** 2
            * c ** 2
            * (
                a ** 2 * b ** 2 * m ** 2
                + a ** 2 * c ** 2 * k ** 2
                - a ** 2 * k ** 2 * z ** 2
                + 2 * a ** 2 * k * m * y * z
                - a ** 2 * m ** 2 * y ** 2
                + b ** 2 * c ** 2
                - b ** 2 * m ** 2 * x ** 2
                + 2 * b ** 2 * m * x * z
                - b ** 2 * z ** 2
                - c ** 2 * k ** 2 * x ** 2
                + 2 * c ** 2 * k * x * y
                - c ** 2 * y ** 2
            )
        )
    ) / (a ** 2 * b ** 2 * m ** 2 + a ** 2 * c ** 2 * k ** 2 + b ** 2 * c ** 2)
    p0, p1 = [x + t0, y + k * t0, z + m * t0], [x - t1, y - t1 * k, z - t1 * m]

    return p0, p1


def project_point_on_ellipsoid(x, y, z, a, b, c):
    """
    Return the projection of a point on an ellipsoid.

    Parameters
    ----------
    x, y, z: Cartesian coordinates of point

    a, b, c: semi-axises of the ellipsoid
    """
    p1, p2 = intersection_ellipsoid_line(x, y, z, x, y, z, a, b, c)

    if np.linalg.norm([p1[0] - x, p1[1] - y, p1[2] - z]) <= np.linalg.norm(
        [p2[0] - x, p2[1] - y, p2[2] - z]
    ):
        return p1
    else:
        return p2
