"""Function helpers.

"""
import numpy as np
from astropy import units as u
from astropy.time import Time

from poliastro.core.util import (
    circular_velocity as circular_velocity_fast,
    norm as norm_fast,
    rotate as rotate_fast,
)

u.kms = u.km / u.s
u.km3s2 = u.km ** 3 / u.s ** 2


def circular_velocity(k, a):
    """Compute circular velocity for a given body (k) and semimajor axis (a).

    """
    return circular_velocity_fast(k.to(u.km3s2).value, a.to(u.km).value) * u.kms


def rotate(vector, angle, axis="z"):
    """Rotates a vector around axis a right-handed positive angle.

    This is just a convenience function around
    :py:func:`astropy.coordinates.matrix_utilities.rotation_matrix`.

    Parameters
    ----------
    vector : ~astropy.units.Quantity
        Dimension 3 vector.
    angle : ~astropy.units.Quantity
        Angle of rotation.
    axis : str, optional
        Either 'x', 'y' or 'z'.

    Note
    -----
    This performs a so-called *active* or *alibi* transformation: rotates the
    vector while the coordinate system remains unchanged. To do the opposite
    operation (*passive* or *alias* transformation) call the function as
    ``rotate(vec, ax, -angle, unit)`` or use the convenience function
    :py:func:`transform`, see [1]_.

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Rotation_matrix#Ambiguities

    """
    return (
        rotate_fast(vector.value, angle.to(u.rad).value, ["x", "y", "z"].index(axis))
        * vector.unit
    )


def transform(vector, angle, axis="z"):
    """Rotates a coordinate system around axis a positive right-handed angle.

    Note
    -----
    This is a convenience function, equivalent to
    ``rotate(vec, -angle, axis, unit)``.
    Refer to the documentation of :py:func:`rotate` for further information.

    """
    return rotate(vector, -angle, axis)


def norm(vec):
    """Norm of a Quantity vector that respects units.

    Parameters
    ----------
    vec : ~astropy.units.Quantity
        Vector with units.

    """
    return norm_fast(vec.value) * vec.unit


def time_range(start, *, periods=50, spacing=None, end=None, format=None, scale=None):
    """Generates range of astronomical times.

    .. versionadded:: 0.8.0

    Parameters
    ----------
    periods : int, optional
        Number of periods, default to 50.
    spacing : Time or Quantity, optional
        Spacing between periods, optional.
    end : Time or equivalent, optional
        End date.

    Returns
    -------
    Time
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


@u.quantity_input(ecc=u.one)
def hyp_nu_limit(ecc, r_max_ratio=np.inf):
    r"""Limit true anomaly for hyperbolic orbits.

    Parameters
    ----------
    ecc : ~astropy.units.Quantity
        Eccentricity, should be larger than 1.
    r_max_ratio : float, optional
        Value of :math:`r_{\text{max}} / p` for this angle, default to infinity.

    """
    return np.arccos(-(1 - 1 / r_max_ratio) / ecc)
