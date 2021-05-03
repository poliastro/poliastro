"""Function helpers.

"""
import numpy as np
from astropy import units as u
from astropy.time import Time
from numpy.linalg import norm as norm_np

from .core.util import alinspace as alinspace_fast


def norm(vec):
    """Norm of a Quantity vector that respects units.

    Parameters
    ----------
    vec : ~astropy.units.Quantity
        Vector with units.

    """
    return norm_np(vec.value) * vec.unit


def time_range(start, *, periods=50, spacing=None, end=None, format=None, scale=None):
    """Generates range of astronomical times.

    .. versionadded:: 0.8.0

    Parameters
    ----------
    periods : int, optional
        Number of periods, default to 50.
    spacing : ~astropy.time.Time or ~astropy.units.Quantity, optional
        Spacing between periods, optional.
    end : ~astropy.time.Time or equivalent, optional
        End date.

    Returns
    -------
    result: ~astropy.time.Time
        Array of time values.

    """
    start = Time(start, format=format, scale=scale)

    if spacing is not None and end is None:
        result = start + spacing * np.arange(0, periods)

    elif end is not None and spacing is None:
        end = Time(end, format=format, scale=scale)
        result = start + (end - start) * np.linspace(0, 1, periods)

    else:
        raise ValueError("Either 'end' or 'spacing' must be specified")

    return result


@u.quantity_input(value=u.rad, values=u.rad)
def find_closest_value(value, values):
    """Calculates the closest value in the given values.

    Parameters
    ----------
    value : ~astropy.units.Quantity
        Reference value.
    values : ~astropy.units.Quantity
        Values to search from.

    """
    index = np.abs(np.asarray(values) * u.rad - value).argmin()
    return values[index]


@u.quantity_input(start=u.rad, stop=u.rad)
def alinspace(start, stop=None, *, num=50, endpoint=True):
    """Return increasing, evenly spaced angular values over a specified interval."""
    if stop is None:
        stop = start + 2 * np.pi * u.rad

    return (
        alinspace_fast(start.to_value(u.rad), stop.to_value(u.rad), num, endpoint)
        * u.rad
    )
