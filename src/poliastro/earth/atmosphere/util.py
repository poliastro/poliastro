""" This script holds several utilities related to atmospheric computations. """


def geometric_to_geopotential(z, r0):
    """Converts from given geometric altitude to geopotential one.

    Parameters
    ----------
    z: ~astropy.units.Quantity
        Geometric altitude.
    r0: ~astropy.units.Quantity
        Planet/Natural satellite radius.

    Returns
    -------
    h: ~astropy.units.Quantity
        Geopotential altitude.
    """

    h = r0 * z / (r0 + z)
    return h


z_to_h = geometric_to_geopotential


def geopotential_to_geometric(h, r0):
    """Converts from given geopotential altitude to geometric one.

    Parameters
    ----------
    h: ~astropy.units.Quanity
        Geopotential altitude.
    r0: ~astropy.units.Quantity
        Planet/Natural satellite radius.

    Returns
    -------
    z: ~astropy.units.Quantity
        Geometric altitude.
    """

    z = r0 * h / (r0 - h)
    return z


h_to_z = geopotential_to_geometric


def gravity(z, g0, r0):
    """Relates Earth gravity field magnitude with the geometric height.

    Parameters
    ----------
    z: ~astropy.units.Quantity
        Geometric height.
    g0: ~astropy.units.Quantity
        Gravity value at sea level.
    r0: ~astropy.units.Quantity
        Planet/Natural satellite radius.

    Returns
    -------
    g: ~astropy.units.Quantity
        Gravity value at given geometric altitude.
    """

    g = g0 * (r0 / (r0 + z)) ** 2
    return g
