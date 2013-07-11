"""Function helpers.

"""

import numpy as np
from numpy import cos, sin, pi

__all__ = ['rotate']

twopi = 2 * pi


def rotate(vec, ax, angle):
    """Rotates a vector around axis a right-handed positive angle.

    This performs a so-called active or alibi transformation: rotates the
    vector while the coordinate system remains unchanged. To do the opposite
    operation (passive or alias transformation) call the function as
    ``rotate(vec, ax, -angle)`` or use the convenience function `transform`.

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
    Based on
    http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle

    """
    # TODO: Broadcast angle and axis too, will have to do it OVER THE THIRD DIMENSION LOL
    #assert ax.shape == (3,)
    #assert vec.shape[-1] == 3
    if type(ax) is int:
        ax_ = np.zeros(3)
        ax_[ax - 1] = 1
        ax = ax_
    rot = np.zeros((3, 3))
    rot[..., :, :] = np.array([
        [
            cos(angle) + ax[0] ** 2 * (1 - cos(angle)),
            ax[0] * ax[1] * (1 - cos(angle)) - ax[2] * sin(angle),
            ax[0] * ax[2] * (1 - cos(angle)) + ax[1] * sin(angle)
        ],
        [
            ax[1] * ax[0] * (1 - cos(angle)) + ax[2] * sin(angle),
            cos(angle) + ax[1] ** 2 * (1 - cos(angle)),
            ax[1] * ax[2] * (1 - cos(angle)) - ax[0] * sin(angle)
        ],
        [
            ax[2] * ax[0] * (1 - cos(angle)) - ax[1] * sin(angle),
            ax[2] * ax[1] * (1 - cos(angle)) + ax[0] * sin(angle),
            cos(angle) + ax[2] ** 2 * (1 - cos(angle))]
    ])
    return np.dot(rot, vec)


def transform(vec, ax, angle):
    """Rotates a coordinate system around axis a positive right-handed angle.

    Notes
    -----
    This is a convenience function, equivalent to `rotate(vec, ax, -angle)`.
    Refer to the documentation of that function for further information.

    """
    return rotate(vec, ax, -angle)


def normalize(angle):
    """Normalize angle between 0 and 2 pi radians.

    """
    return angle % twopi
