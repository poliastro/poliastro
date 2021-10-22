import numpy as np
from numba import njit as jit
from numpy import cos, sin


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
