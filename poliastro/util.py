"""Function helpers.

"""

import numpy as np
from numpy import cos, sin, pi

__all__ = ['rotate', 'transform', 'normalize', 'direct_angles']


def rotate(vec, ax, angle):
    """Rotates a vector around axis a right-handed positive angle.

    Parameters
    ----------
    vec : array
        Dimension 3 vector.
    ax : array, int
        Axis of rotation. If 1, 2, 3 is given, then one of the coordinate axis
        is chosen.
    angle : float
        Angle of rotation (rad).

    Notes
    -----
    This performs a so-called active or alibi transformation: rotates the
    vector while the coordinate system remains unchanged. To do the opposite
    operation (passive or alias transformation) call the function as
    ``rotate(vec, ax, -angle)`` or use the convenience function `transform`,
    see `[1]_`.

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Rotation_matrix#Ambiguities

    """
    if type(ax) is int:
        ax_ = np.zeros(3)
        ax_[ax - 1] = 1
        ax = ax_
    rot = _rot_matrix(ax, angle)
    return np.dot(rot, vec)


def transform(vec, ax, angle):
    """Rotates a coordinate system around axis a positive right-handed angle.

    Notes
    -----
    This is a convenience function, equivalent to `rotate(vec, ax, -angle)`.
    Refer to the documentation of that function for further information.

    """
    return rotate(vec, ax, -angle)


def _rot_matrix(ax, ang):
    """
    Based on
    http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle

    """
    rot = np.array([
        [
            cos(ang) + ax[0] ** 2 * (1 - cos(ang)),
            ax[0] * ax[1] * (1 - cos(ang)) - ax[2] * sin(ang),
            ax[0] * ax[2] * (1 - cos(ang)) + ax[1] * sin(ang)
            ],
        [
            ax[1] * ax[0] * (1 - cos(ang)) + ax[2] * sin(ang),
            cos(ang) + ax[1] ** 2 * (1 - cos(ang)),
            ax[1] * ax[2] * (1 - cos(ang)) - ax[0] * sin(ang)
            ],
        [
            ax[2] * ax[0] * (1 - cos(ang)) - ax[1] * sin(ang),
            ax[2] * ax[1] * (1 - cos(ang)) + ax[0] * sin(ang),
            cos(ang) + ax[2] ** 2 * (1 - cos(ang))]
    ])
    return rot


def normalize(angle):
    """Normalize angle between 0 and 2 pi radians.

    """
    twopi = 2 * np.pi
    return angle % twopi


def direct_angles(a_i, a_f):
    """Return angles so a_f > a_i, with a1 unchanged.

    """
    twopi = 2 * np.pi
    a_f = a_f - twopi * ((a_f - a_i) // twopi)
    return a_i, a_f
