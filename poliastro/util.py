"""Function helpers.

"""

import numpy as np
from numpy import cos, sin

__all__ = ['rotate']


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
