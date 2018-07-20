import numpy as np
from numpy import sin, cos

from poliastro.core.jit import jit

@jit
def rotate(vec, ax, angle):
    """Rotates the coordinate system around axis x, y or z a CCW angle.
    Parameters
    ----------
    vec : ndarray
        Dimension 3 vector.
    ax : int
        Axis to be rotated.
    angle : float
        Angle of rotation (rad).
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
    if ax == 0:
        sl = slice(1, 3, 1)
    elif ax == 1:
        sl = slice(0, 3, 2)
    elif ax == 2:
        sl = slice(0, 2, 1)
    else:
        raise ValueError("Invalid axis: must be one of 0, 1 or 2")

    rot[sl, sl] = np.array((
        (cos(angle), -sin(angle)),
        (sin(angle), cos(angle))
    ))
    if ax == 1:
        rot = rot.T

    # np.dot() arguments must all have the same dtype
    return np.dot(rot, vec.astype(rot.dtype))


@jit
def transform(vec, ax, angle):
    """Rotates a coordinate system around axis a positive right-handed angle.
    Notes
    -----
    This is a convenience function, equivalent to `rotate(vec, ax, -angle)`.
    Refer to the documentation of that function for further information.
    """
    return rotate(vec, ax, -angle)
