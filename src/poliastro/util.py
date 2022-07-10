"""Function helpers.

"""
import numpy as np
from astropy import units as u
from astropy.coordinates import Angle
from astropy.time import Time

from poliastro._math.linalg import norm as norm_fast
from poliastro.core.util import alinspace as alinspace_fast


def norm(vec, axis=None):
    """Norm of a Quantity vector that respects units.

    Parameters
    ----------
    vec : ~astropy.units.Quantity
        Vector with units.
    axis : int or None
        Axis along which to compute the vector norms.

    """
    if axis is not None:
        # axis keyword argument is not yet supported by numba,
        # see https://github.com/numba/numba/issues/2558
        from numpy.linalg import norm as norm_np

        result = norm_np(vec.value, axis=axis)
    else:
        result = norm_fast(vec.value)

    return result << vec.unit


def time_range(
    start, *, periods=50, spacing=None, end=None, format=None, scale=None
):
    """Generates range of astronomical times.

    .. versionadded:: 0.8.0

    Parameters
    ----------
    start : ~astropy.time.Time or ~astropy.units.Quantity
        Start time.
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
        alinspace_fast(
            start.to_value(u.rad), stop.to_value(u.rad), num, endpoint
        )
        * u.rad
    )


@u.quantity_input(angle=u.rad, limit=u.rad)
def wrap_angle(angle, limit=180 * u.deg):
    return Angle(angle).wrap_at(limit)
