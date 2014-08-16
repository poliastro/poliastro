# coding: utf-8
"""Function helpers.

"""

import numpy as np

from astropy.coordinates import angles
from astropy import units as u
u.one = u.dimensionless_unscaled  # astropy #1980


def rotate(vector, angle, axis='z', unit=None):
    """Rotates a vector around axis a right-handed positive angle.

    This is just a convenience function around
    `astropy.coordinates.angles.rotation_matrix`.

    Parameters
    ----------
    vector : array_like
        Dimension 3 vector.
    angle : convertible to Angle
        Angle of rotation.
    axis : str or 3-sequence
        Either 'x','y', 'z', or a (x,y,z) specifying an axis to rotate
        about. If 'x','y', or 'z', the rotation sense is
        counterclockwise looking down the + axis (e.g. positive
        rotations obey left-hand-rule).
    unit : UnitBase, optional
        If `angle` does not have associated units, they are in this
        unit.  If neither are provided, it is assumed to be degrees.

    Notes
    -----
    This is just a convenience function around
    `astropy.coordinates.angles.rotation_matrix`.
    This performs a so-called *active* or *alibi* transformation: rotates the
    vector while the coordinate system remains unchanged. To do the opposite
    operation (*passive* or *alias* transformation) call the function as
    ``rotate(vec, ax, -angle, unit)`` or use the convenience function
    :py:func:`transform`, see [1].

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Rotation_matrix#Ambiguities

    """
    rot = np.asarray(angles.rotation_matrix(-angle, axis, unit))
    return np.dot(rot, vector)


def transform(vector, angle, axis='z', unit=None):
    """Rotates a coordinate system around axis a positive right-handed angle.

    Notes
    -----
    This is a convenience function, equivalent to
    ``rotate(vec, -angle, axis, unit)``.
    Refer to the documentation of :py:func:`rotate` for further information.

    """
    return rotate(vector, -angle, axis, unit)


def check_units(quantities, units):
    """Check if list of quantities is consistent with given list of units.

    """
    return all(unit.is_equivalent((quant * u.one).unit)
               for unit, quant in zip(units, quantities))


def norm(vec):
    """Norm of a Quantity vector that respects units.

    """
    return np.sqrt(vec.dot(vec))
