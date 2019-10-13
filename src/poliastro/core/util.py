import numpy as np
from numpy import cos, sin

from ._jit import jit


@jit
def circular_velocity(k, a):
    r"""Compute circular velocity for a given body given thegravitational parameter and the semimajor axis.

    .. math::

       v = \sqrt{\frac{\mu}{a}}

    Parameters
    ----------

    k : float
        Gravitational Parameter
    a : float
        Semimajor Axis

    """
    return np.sqrt(k / a)


@jit
def rotate(vec, angle, axis):
    r"""Rotates a vector around axis x, y or z a Counter ClockWise angle.

    .. math::

        R(\theta) = \begin{bmatrix}
            \cos(\theta) & -\sin(\theta) \\
            \sin(\theta) & \cos(\theta)
            \end{bmatrix}

    Parameters
    ----------
    vec : ndarray
        Dimension 3 vector.
    angle : float
        Angle of rotation (rad).
    axis : int
        Axis to be rotated (X-axis: 0, Y-axis: 1, Z-axis: 2)

    Note
    ----
    This performs a so-called active or alibi transformation: rotates the
    vector while the coordinate system remains unchanged. To do the opposite
    operation (passive or alias transformation) call the function as
    `rotate(vec, ax, -angle)` or use the convenience function `transform`.

    References
    ----------
    http://en.wikipedia.org/wiki/Rotation_matrix#Ambiguities

    """

    assert vec.shape == (3,)

    rot = np.eye(3)
    if axis == 0:
        sl = slice(1, 3, 1)
    elif axis == 1:
        sl = slice(0, 3, 2)
        angle = -angle
    elif axis == 2:
        sl = slice(0, 2, 1)
    else:
        raise ValueError("Invalid axis: must be one of 'x', 'y' or 'z'")
    rot[sl, sl] = np.array(((cos(angle), -sin(angle)), (sin(angle), cos(angle))))

    # np.dot() arguments must all have the same dtype
    return np.dot(rot, vec.astype(rot.dtype))


@jit
def transform(vec, angle, axis):
    """Rotates a coordinate system around axis a positive right-handed angle.

    Note
    -----
    This is a convenience function, equivalent to `rotate(vec, ax, -angle)`.
    Refer to the documentation of that function for further information.

    """
    return rotate(vec, -angle, axis)


@jit
def norm(vec):
    r"""Returns the norm of a 3 dimension vector.

    .. math::

        \left \| \vec{v} \right \| = \sqrt{\sum_{i=1}^{n}v_{i}^2}

    Parameters
    ----------
    vec: ndarray
        Dimension 3 vector.

    Examples
    --------
    >>> vec = np.array([1, 1, 1])
    >>> norm(vec)
    1.7320508075688772

    """
    vec = 1.0 * vec  # Cast to float
    return np.sqrt(vec.dot(vec))


@jit
def rotation_matrix(angle, axis):
    c = cos(angle)
    s = sin(angle)
    if axis == 0:
        return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])
    elif axis == 1:
        return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [s, 0.0, c]])
    elif axis == 2:
        return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    else:
        raise ValueError("Invalid axis: must be one of 'x', 'y' or 'z'")
