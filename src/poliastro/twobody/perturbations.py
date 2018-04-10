from astropy import units as u
import numpy as np
from poliastro.util import norm
from poliastro.bodies import Sun, Earth


def J2_perturbation(t0, state, k, J2, R):
    """Calculates J2_perturbation acceleration (km/s2)

    Parameters
    ----------
    t0 : float
        Current time (s)
    state : numpy.ndarray
        Six component state vector [x, y, z, vx, vy, vz] (km, km/s).
    k : float
        gravitational constant, (km^3/s^2)
    J2: float
        obliteness factor
    R: float
        attractor radius

    Notes
    -----
    The J2 accounts for the obliteness of the attractor. The formula is given in
    Howard Curtis, (12.30)
    """

    r_vec = state[:3]
    r = norm(r_vec)

    factor = (3.0 / 2.0) * k * J2 * (R ** 2) / (r ** 5)

    a_x = 5.0 * r_vec[2] ** 2 / r ** 2 - 1
    a_y = 5.0 * r_vec[2] ** 2 / r ** 2 - 1
    a_z = 5.0 * r_vec[2] ** 2 / r ** 2 - 3
    return np.array([a_x, a_y, a_z]) * r_vec * factor


def atmospheric_drag(t0, state, k, R, B, H0, rho0):
    """Calculates atmospheric drag acceleration (km/s2)

    Parameters
    ----------
    t0 : float
        Current time (s)
    state : numpy.ndarray
        Six component state vector [x, y, z, vx, vy, vz] (km, km/s).
    k : float
        gravitational constant, (km^3/s^2)
    B: float
        ballistic coefficient B = C_D[] x A[m^2] / m[kg], (m^2 / kg)
    H0 : float
        the e-times density decay heigth, (km)
    rho0: float
        the exponent density pre-factor, (kg / m^3)

    Notes
    -----
    This function provides the acceleration due to atmospheric drag. We follow
    Howard Curtis, section 12.4
    the atmospheric density model is rho(H) = rho0 x exp(-H / H0)
    """

    H = norm(state[:3])

    v_vec = state[3:]
    v = norm(v_vec)

    rho = rho0 * np.exp(-(H - R) / H0)

    return -(1.0 / 2.0) * rho * B * v * v_vec
