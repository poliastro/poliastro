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


def intersection_ellipsoid_line(a, b, c, u1, u2, u3, x, y, z):
    """
    Intersection of an ellipsoid defined by its axises a, b, c with the
    line p + λu

    Parameters
    ----------
    a, b, c: float
        the ellipsoidal axises
    u1, u2, u3: float
        the line vector
    x, y, z: float
        a point of the line

    Returns
    -------
    This returns both of the points intersecting the ellipsoid.

    """
    # Get rid of one parameter by translating the line's direction vector
    k, l = u2 / u1, u3 / u1
    t0 = (
        -np.power(a, 2) * np.power(b, 2) * l * z
        - np.power(a, 2) * np.power(c, 2) * k * y
        - np.power(b, 2) * np.power(c, 2) * x
        + np.sqrt(
            np.power(a, 2)
            * np.power(b, 2)
            * np.power(c, 2)
            * (
                np.power(a, 2) * np.power(b, 2) * np.power(l, 2)
                + np.power(a, 2) * np.power(c, 2) * np.power(k, 2)
                - np.power(a, 2) * np.power(k, 2) * np.power(z, 2)
                + 2 * np.power(a, 2) * k * l * y * z
                - np.power(a, 2) * np.power(l, 2) * np.power(y, 2)
                + np.power(b, 2) * np.power(c, 2)
                - np.power(b, 2) * np.power(l, 2) * np.power(x, 2)
                + 2 * np.power(b, 2) * l * x * z
                - np.power(b, 2) * np.power(z, 2)
                - np.power(c, 2) * np.power(k, 2) * np.power(x, 2)
                + 2 * np.power(c, 2) * k * x * y
                - np.power(c, 2) * np.power(y, 2)
            )
        )
    ) / (
        np.power(a, 2) * np.power(b, 2) * np.power(l, 2)
        + np.power(a, 2) * np.power(c, 2) * np.power(k, 2)
        + np.power(b, 2) * np.power(c, 2)
    )
    t1 = (
        np.power(a, 2) * np.power(b, 2) * l * z
        + np.power(a, 2) * np.power(c, 2) * k * y
        + np.power(b, 2) * np.power(c, 2) * x
        + np.sqrt(
            np.power(a, 2)
            * np.power(b, 2)
            * np.power(c, 2)
            * (
                np.power(a, 2) * np.power(b, 2) * np.power(l, 2)
                + np.power(a, 2) * np.power(c, 2) * np.power(k, 2)
                - np.power(a, 2) * np.power(k, 2) * np.power(z, 2)
                + 2 * np.power(a, 2) * k * l * y * z
                - np.power(a, 2) * np.power(l, 2) * np.power(y, 2)
                + np.power(b, 2) * np.power(c, 2)
                - np.power(b, 2) * np.power(l, 2) * np.power(x, 2)
                + 2 * np.power(b, 2) * l * x * z
                - np.power(b, 2) * np.power(z, 2)
                - np.power(c, 2) * np.power(k, 2) * np.power(x, 2)
                + 2 * np.power(c, 2) * k * x * y
                - np.power(c, 2) * np.power(y, 2)
            )
        )
    ) / (
        np.power(a, 2) * np.power(b, 2) * np.power(l, 2)
        + np.power(a, 2) * np.power(c, 2) * np.power(k, 2)
        + np.power(b, 2) * np.power(c, 2)
    )
    p0, p1 = [x + t0, y + k * t0, z + l * t0], [x - t1, y - t1 * k, z - t1 * l]

    return p0, p1


def magnitude(u, v):
    return np.sqrt(sum([(a - b) ** 2 for a, b in zip(u, v)]))
