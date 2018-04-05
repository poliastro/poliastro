from astropy import units as u
import numpy as np
from poliastro.util import norm
from poliastro.bodies import Sun, Earth


def J2_perturbation(t0, state, k, J2, R):
    r"""Calculates J2_perturbation acceleration (km/s2)

    Parameters
    ----------
    t0 : float
        Current time (s)
    u : numpy.ndarray
        Six component state vector [x, y, z, vx, vy, vz] (km, km/s).
    k : float
        Gravitational constant of main attractor (km^3 / s^2),
    body : float, optional
        The attracting body, default to Earth

    Notes
    -----
    The J2 accounts for the obliteness of the attractor. The formula is given in
    Howard Curtis, (12.30)

    """

    r_vec = state[:3]
    r = norm(r_vec)

    factor = (3.0 / 2.0) * k * (R.to(u.km).value ** 2) * J2 / (r ** 5)

    a_x = 5.0 * r_vec[2] ** 2 / r ** 2 - 1
    a_y = 5.0 * r_vec[2] ** 2 / r ** 2 - 1
    a_z = 5.0 * r_vec[2] ** 2 / r ** 2 - 3
    return np.array([a_x, a_y, a_z]) * r_vec * factor
