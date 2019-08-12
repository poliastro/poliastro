""" This script holds several utilities related to atmospheric computations. """
from astropy import units as u
from astropy.constants import R as R_u, g0

# ATMOSPHERIC AND AIR PROPETIES
M0_air = 28.964420 * u.kg / u.kmol
R_air = R_u / M0_air


def geometric_to_geopotential(h, R_body):
    """ Converts from given geometric altitude to geopotential one.

    Parameters
    ----------
    h: ~astropy.units.Quantity
        Geometric altitude.
    R_body: ~astropy.units.Quantity
        Planet/Natural satellite radius.

    Returns
    -------
    H: ~astropy.units.Quantity
        Geopotential altitude.
    """

    H = R_body * h / (R_body + h)
    return H


def geopotential_to_geometric(H, R_body):
    """ Converts from given geopotential altitude to geometric one.

    Parameters
    ----------
    H: ~astropy.units.Quanity
        Geopotential altitude.
    R_body: ~astropy.units.Quantity
        Planet/Natural satellite radius.

    Returns
    -------
    h: ~astropy.units.Quantity
        Geometric altitude.
    """

    h = R_body * H / (R_body - H)
    return h


def gravity(Z, R_body):
    """ Relates Earth gravity field magnitude with the geometric height.

    Parameters
    ----------
    Z: ~astropy.units.Quantity
        Geometric height.
    R_body: ~astropy.units.Quantity
        Planet/Natural satellite radius.

    Returns
    -------
    g: ~astropy.units.Quantity
        Earth's gravity field magnitude.
    """

    g = g0 * (R_body / (R_body + Z)) ** 2
    return g
