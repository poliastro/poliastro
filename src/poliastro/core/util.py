"""Function helpers.

"""
import numpy as np
from numpy import sin, cos

from ._jit import jit


@jit
def circular_velocity(k, a):
    """Compute circular velocity for a given body (k) and semimajor axis (a).

    """
    return np.sqrt(k / a)


@jit
def rotate(vec, angle, axis):
    """Rotates the coordinate system around axis x, y or z a CCW angle.

    Parameters
    ----------
    vec : ndarray
        Dimension 3 vector.
    angle : float
        Angle of rotation (rad).
    axis : int
        Axis to be rotated.

    Notes
    -----
    This performs a so-called active or alibi transformation: rotates the
    vector while the coordinate system remains unchanged. To do the opposite
    operation (passive or alias transformation) call the function as
    `rotate(vec, ax, -angle)` or use the convenience function `transform`,
    see `[1]_`.

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Rotation_matrix#Ambiguities

    """
    assert vec.shape == (3,)

    rot = np.eye(3)
    if axis == 0:
        sl = slice(1, 3, 1)
    elif axis == 1:
        sl = slice(0, 3, 2)
    elif axis == 2:
        sl = slice(0, 2, 1)
    else:
        raise ValueError("Invalid axis: must be one of 'x', 'y' or 'z'")

    rot[sl, sl] = np.array((
        (cos(angle), -sin(angle)),
        (sin(angle), cos(angle))
    ))
    if axis == 1:
        rot = rot.T

    # np.dot() arguments must all have the same dtype
    return np.dot(rot, vec.astype(rot.dtype))


@jit
def transform(vec, angle, axis):
    """Rotates a coordinate system around axis a positive right-handed angle.

    Notes
    -----
    This is a convenience function, equivalent to `rotate(vec, ax, -angle)`.
    Refer to the documentation of that function for further information.

    """
    return rotate(vec, -angle, axis)


@jit
def norm(vec):
    """Norm of a 3d vector.

    """
    vec = 1.0 * vec  # Cast to float
    return np.sqrt(vec.dot(vec))


@jit
def cross(a, b):
    """Computes cross product between two vectors"""
    # np.cross is not supported in numba nopython mode, see
    # https://github.com/numba/numba/issues/2978
    return np.array((
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    ))
