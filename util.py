"""Function helpers

"""

import numpy as np

__all__ = ['days_from_epoch', 'rotate']


def days_from_epoch(d):
    # TODO: Docstring
    delta = d - EPOCH
    return delta.days


def rotate(vec, ax, angle):
    """Rotates the coordinate system around axis 1, 2 or 3 a CCW angle.

    Parameters
    ----------
    vec : array
        Dimension 3 vector.
    ax : int
        Axis to be rotated.
    angle : float
        Angle of rotation (rad).

    """
    # TODO: Accept arbitrary axis
    # http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
    assert vec.shape == (3,)
    rot = np.eye(3)
    if ax == 1:
        sl = slice(1, 3)
    elif ax == 2:
        sl = slice(0, 3, 2)
    elif ax == 3:
        sl = slice(0, 2)
    rot[sl, sl] = np.array([
        [np.cos(angle), np.sin(angle)],
        [-np.sin(angle), np.cos(angle)]
    ])
    return np.dot(rot, vec)
