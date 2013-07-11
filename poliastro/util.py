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
    angle : array
        Angle(s) of rotation (rad).

    Notes
    -----
    This function performs some unholy sorcery to work also with an array of
    vectors and an array of angles. There be dragons!

    """
    if type(ax) is int:
        ax_ = np.zeros(3)
        ax_[ax - 1] = 1
        ax = ax_
    angle = np.atleast_1d(angle)
    vec = np.atleast_2d(vec.T).T
    angle, vec = np.broadcast_arrays(angle, vec)
    rot = _rot_matrix(ax, angle[0])
    qqq = np.einsum('aik,km->aim', rot, vec)
    rot_vec = np.einsum('lhl->hl', qqq)
    if rot_vec.shape[-1] == 1:
        rot_vec = rot_vec[:, 0]
    return rot_vec


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
    rot = np.zeros((3, 3) + ang.shape)
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
    return rot.transpose(2, 0, 1)


def normalize(angle):
    """Normalize angle between 0 and 2 pi radians.

    """
    return angle % twopi


def direct_angles(a1, a2):
    twopi = 2 * np.pi
    a2 = a2 - twopi * ((a2 - a1) // twopi)
    return a1, a2
