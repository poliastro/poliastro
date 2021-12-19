import numpy as np
from numba import njit as jit
from numpy import cos, sin
from numpy.linalg import norm as norm


@jit
def eccentricity_vector(k, r, v):
    return ((v @ v - k / norm(r)) * r - (r @ v) * v) / k


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
def rotation_matrix(angle, axis):
    assert axis in (0, 1, 2)
    angle = np.asarray(angle)
    c = cos(angle)
    s = sin(angle)

    a1 = (axis + 1) % 3
    a2 = (axis + 2) % 3
    R = np.zeros(angle.shape + (3, 3))
    R[..., axis, axis] = 1.0
    R[..., a1, a1] = c
    R[..., a1, a2] = -s
    R[..., a2, a1] = s
    R[..., a2, a2] = c
    return R


@jit
def alinspace(start, stop=None, num=50, endpoint=True):
    """Return increasing, evenly spaced angular values over a specified interval."""
    # Create a new variable to avoid numba crash,
    # See https://github.com/numba/numba/issues/5661
    if stop is None:
        stop_ = start + 2 * np.pi
    elif stop <= start:
        stop_ = stop + (np.floor((start - stop) / 2 / np.pi) + 1) * 2 * np.pi
    else:
        stop_ = stop

    if endpoint:
        return np.linspace(start, stop_, num)
    else:
        return np.linspace(start, stop_, num + 1)[:-1]


@jit
def spherical_to_cartesian(v):
    r"""Compute cartesian coordinates from spherical coordinates (norm, colat, long). This function is vectorized.

    .. math::

       v = norm \cdot \begin{bmatrix}
       \sin(colat)\cos(long)\\
       \sin(colat)\sin(long)\\
       \cos(colat)\\
       \end{bmatrix}

    Parameters
    ----------
    v : np.array
        Spherical coordinates in 3D (norm, colat, long). Angles must be in radians.

    Returns
    -------

    v : np.array
        Cartesian coordinates (x,y,z)

    """
    v = np.asarray(v)
    norm = np.expand_dims(np.asarray(v[..., 0]), -1)
    vsin = np.sin(v[..., 1:3])
    vcos = np.cos(v[..., 1:3])
    x = np.asarray(vsin[..., 0] * vcos[..., 1])
    y = np.asarray(vsin[..., 0] * vsin[..., 1])
    z = np.asarray(vcos[..., 0])
    return norm * np.stack((x, y, z), axis=-1)
