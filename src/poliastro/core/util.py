from numba import njit as jit
import numpy as np
from numpy import cos, sin


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
    v : numpy.ndarray
        Spherical coordinates in 3D (norm, colat, long). Angles must be in radians.

    Returns
    -------
    v : numpy.ndarray
        Cartesian coordinates (x,y,z)

    """
    v = np.asarray(v)
    norm_vecs = np.expand_dims(np.asarray(v[..., 0]), -1)
    vsin = np.sin(v[..., 1:3])
    vcos = np.cos(v[..., 1:3])
    x = np.asarray(vsin[..., 0] * vcos[..., 1])
    y = np.asarray(vsin[..., 0] * vsin[..., 1])
    z = np.asarray(vcos[..., 0])
    return norm_vecs * np.stack((x, y, z), axis=-1)


@jit
def planetocentric_to_AltAz(theta, phi):
    r"""Defines transformation matrix to convert from Planetocentric coordinate system
    to the Altitude-Azimuth system.

    .. math::
       t\_matrix = \begin{bmatrix}
       -\sin(theta) & \cos(theta) & 0\\
       -\sin(phi)\cdot\cos(theta) & -\sin(phi)\cdot\sin(theta) & \cos(phi)\\
       \cos(phi)\cdot\cos(theta) & \cos(phi)\cdot\sin(theta) & \sin(phi)
       \end{bmatrix}

    Parameters
    ----------
    theta: float
        Local sidereal time
    phi: float
        Planetodetic latitude

    Returns
    -------
    t_matrix: numpy.ndarray
        Transformation matrix
    """
    # Transformation matrix for converting planetocentric equatorial coordinates to topocentric horizon system.
    t_matrix = np.array(
        [
            [-np.sin(theta), np.cos(theta), 0],
            [
                -np.sin(phi) * np.cos(theta),
                -np.sin(phi) * np.sin(theta),
                np.cos(phi),
            ],
            [
                np.cos(phi) * np.cos(theta),
                np.cos(phi) * np.sin(theta),
                np.sin(phi),
            ],
        ]
    )
    return t_matrix
