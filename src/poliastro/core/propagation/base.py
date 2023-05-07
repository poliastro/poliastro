from numba import njit as jit
import numpy as np


@jit
def func_twobody(t0, u_, k):
    """Differential equation for the initial value two body problem.

    Parameters
    ----------
    t0 : float
        Time.
    u_ : numpy.ndarray
        Six component state vector [x, y, z, vx, vy, vz] (km, km/s).
    k : float
        Standard gravitational parameter.

    """
    x, y, z, vx, vy, vz = u_
    r3 = (x**2 + y**2 + z**2) ** 1.5

    du = np.array([vx, vy, vz, -k * x / r3, -k * y / r3, -k * z / r3])
    return du
