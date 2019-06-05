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
