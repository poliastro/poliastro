""" This script holds several utilities related to atmospheric computations."""

from numba import njit as jit


@jit
def geometric_to_geopotential(z, r0):
    """Converts from given geometric altitude to geopotential one.

    Parameters
    ----------
    z : float
        Geometric altitude.
    r0 : float
        Planet/Natural satellite radius.

    Returns
    -------
    h: float
        Geopotential altitude.
    """

    h = r0 * z / (r0 + z)
    return h


z_to_h = geometric_to_geopotential


@jit
def geopotential_to_geometric(h, r0):
    """Converts from given geopotential altitude to geometric one.

    Parameters
    ----------
    h : float
        Geopotential altitude.
    r0 : float
        Planet/Natural satellite radius.

    Returns
    -------
    z: float
        Geometric altitude.
    """

    z = r0 * h / (r0 - h)
    return z


h_to_z = geopotential_to_geometric


@jit
def gravity(z, g0, r0):
    """Relates Earth gravity field magnitude with the geometric height.

    Parameters
    ----------
    z : float
        Geometric height.
    g0 : float
        Gravity value at sea level.
    r0 : float
        Planet/Natural satellite radius.

    Returns
    -------
    g: float
        Gravity value at given geometric altitude.
    """

    g = g0 * (r0 / (r0 + z)) ** 2
    return g


@jit
def _get_index(x, x_levels):
    """Finds element in list and returns index.

    Parameters
    ----------
    x : float
        Element to be searched.
    x_levels : list
        List for searching.

    Returns
    -------
    i: int
        Index for the value.

    """
    for i, value in enumerate(x_levels):
        if i < len(x_levels) and value < x:
            continue
        elif x == value:
            return i
        else:
            return i - 1


@jit
def _check_altitude(alt, r0, geometric):
    # Get geometric and geopotential altitudes
    if geometric:
        z = alt
        h = z_to_h(z, r0)
    else:
        h = alt
        z = h_to_z(h, r0)

    return z, h
