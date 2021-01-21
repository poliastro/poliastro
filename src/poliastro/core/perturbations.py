import numpy as np
from astropy import units as u
from numba import njit as jit
from numpy.linalg import norm


@jit
def J2_perturbation(t0, state, k, J2, R):
    r"""Calculates J2_perturbation acceleration (km/s2)

    .. math::

        \vec{p} = \frac{3}{2}\frac{J_{2}\mu R^{2}}{r^{4}}\left [\frac{x}{r}\left ( 5\frac{z^{2}}{r^{2}}-1 \right )\vec{i} + \frac{y}{r}\left ( 5\frac{z^{2}}{r^{2}}-1 \right )\vec{j} + \frac{z}{r}\left ( 5\frac{z^{2}}{r^{2}}-3 \right )\vec{k}\right]

    .. versionadded:: 0.9.0

    Parameters
    ----------
    t0 : float
        Current time (s)
    state : numpy.ndarray
        Six component state vector [x, y, z, vx, vy, vz] (km, km/s).
    k : float
        gravitational constant, (km^3/s^2)
    J2: float
        oblateness factor
    R: float
        attractor radius

    Note
    ----
    The J2 accounts for the oblateness of the attractor. The formula is given in
    Howard Curtis, (12.30)

    """
    r_vec = state[:3]
    r = norm(r_vec)

    factor = (3.0 / 2.0) * k * J2 * (R ** 2) / (r ** 5)

    a_x = 5.0 * r_vec[2] ** 2 / r ** 2 - 1
    a_y = 5.0 * r_vec[2] ** 2 / r ** 2 - 1
    a_z = 5.0 * r_vec[2] ** 2 / r ** 2 - 3
    return np.array([a_x, a_y, a_z]) * r_vec * factor


@jit
def J3_perturbation(t0, state, k, J3, R):
    r"""Calculates J3_perturbation acceleration (km/s2)

    Parameters
    ----------
    t0 : float
        Current time (s)
    state : numpy.ndarray
        Six component state vector [x, y, z, vx, vy, vz] (km, km/s).
    k : float
        gravitational constant, (km^3/s^2)
    J3: float
        oblateness factor
    R: float
        attractor radius

    Note
    ----
    The J3 accounts for the oblateness of the attractor. The formula is given in
    Howard Curtis, problem 12.8
    This perturbation has not been fully validated, see https://github.com/poliastro/poliastro/pull/398

    """
    r_vec = state[:3]
    r = norm(r_vec)

    factor = (1.0 / 2.0) * k * J3 * (R ** 3) / (r ** 5)
    cos_phi = r_vec[2] / r

    a_x = 5.0 * r_vec[0] / r * (7.0 * cos_phi ** 3 - 3.0 * cos_phi)
    a_y = 5.0 * r_vec[1] / r * (7.0 * cos_phi ** 3 - 3.0 * cos_phi)
    a_z = 3.0 * (35.0 / 3.0 * cos_phi ** 4 - 10.0 * cos_phi ** 2 + 1)
    return np.array([a_x, a_y, a_z]) * factor


@jit
def atmospheric_drag_exponential(t0, state, k, R, C_D, A_over_m, H0, rho0):
    r"""Calculates atmospheric drag acceleration (km/s2)

    .. math::

        \vec{p} = -\frac{1}{2}\rho v_{rel}\left ( \frac{C_{d}A}{m} \right )\vec{v_{rel}}


    .. versionadded:: 0.9.0

    Parameters
    ----------
    t0 : float
        Current time (s)
    state : numpy.ndarray
        Six component state vector [x, y, z, vx, vy, vz] (km, km/s).
    k : float
        gravitational constant, (km^3/s^2)
    R : float
        radius of the attractor (km)
    C_D: float
        dimensionless drag coefficient ()
    A_over_m: float
        frontal area/mass of the spacecraft (km^2/kg)
    H0 : float
        atmospheric scale height, (km)
    rho0: float
        the exponent density pre-factor, (kg / km^3)

    Note
    ----
    This function provides the acceleration due to atmospheric drag
    using an overly-simplistic exponential atmosphere model. We follow
    Howard Curtis, section 12.4
    the atmospheric density model is rho(H) = rho0 x exp(-H / H0)

    """
    H = norm(state[:3])

    v_vec = state[3:]
    v = norm(v_vec)
    B = C_D * A_over_m
    rho = rho0 * np.exp(-(H - R) / H0)

    return -(1.0 / 2.0) * rho * B * v * v_vec


def atmospheric_drag_model(t0, state, k, R, C_D, A_over_m, model):
    r"""Calculates atmospheric drag acceleration (km/s2)

    .. math::

        \vec{p} = -\frac{1}{2}\rho v_{rel}\left ( \frac{C_{d}A}{m} \right )\vec{v_{rel}}


    .. versionadded:: 1.14

    Parameters
    ----------
    t0 : float
        Current time (s)
    state : numpy.ndarray
        Six component state vector [x, y, z, vx, vy, vz] (km, km/s).
    k : float
        gravitational constant, (km^3/s^2)
    R : float
        radius of the attractor (km)
    C_D: float
        dimensionless drag coefficient ()
    A_over_m: float
        frontal area/mass of the spacecraft (km^2/kg)
    model: a callable model from poliastro.earth.atmosphere

    Note
    ----
    This function provides the acceleration due to atmospheric drag, as
    computed by a model from poliastro.earth.atmosphere

    """
    H = norm(state[:3])

    v_vec = state[3:]
    v = norm(v_vec)
    B = C_D * A_over_m

    if H < R:
        # The model doesn't want to see a negative altitude
        # The integration will go a little negative searching for H = R
        H = R

    rho = model.density((H - R) * u.km).to(u.kg / u.km ** 3).value

    return -(1.0 / 2.0) * rho * B * v * v_vec


@jit
def shadow_function(r_sat, r_sun, R):
    r"""Determines whether the satellite is in attractor's shadow, uses algorithm 12.3 from Howard Curtis

    Parameters
    ----------
    r_sat : numpy.ndarray
        position of the satellite in the frame of attractor (km)
    r_sun : numpy.ndarray
        position of star in the frame of attractor (km)
    R : float
        radius of body (attractor) that creates shadow (km)

    """

    r_sat_norm = np.sqrt(np.sum(r_sat ** 2))
    r_sun_norm = np.sqrt(np.sum(r_sun ** 2))

    theta = np.arccos(np.dot(r_sat, r_sun) / r_sat_norm / r_sun_norm)
    theta_1 = np.arccos(R / r_sat_norm)
    theta_2 = np.arccos(R / r_sun_norm)

    return theta < theta_1 + theta_2


def third_body(t0, state, k, k_third, perturbation_body):
    r"""Calculates 3rd body acceleration (km/s2)

    .. math::

        \vec{p} = \mu_{m}\left ( \frac{\vec{r_{m/s}}}{r_{m/s}^3} - \frac{\vec{r_{m}}}{r_{m}^3} \right )

    Parameters
    ----------
    t0 : float
        Current time (s)
    state : numpy.ndarray
        Six component state vector [x, y, z, vx, vy, vz] (km, km/s).
    k : float
        gravitational constant, (km^3/s^2)
    perturbation_body: a callable object returning the position of the pertubation body that causes the perturbation

    Note
    ----
    This formula is taken from Howard Curtis, section 12.10. As an example, a third body could be
    the gravity from the Moon acting on a small satellite.

    """

    body_r = perturbation_body(t0)
    delta_r = body_r - state[:3]
    return k_third * delta_r / norm(delta_r) ** 3 - k_third * body_r / norm(body_r) ** 3


def radiation_pressure(t0, state, k, R, C_R, A_over_m, Wdivc_s, star):
    r"""Calculates radiation pressure acceleration (km/s2)

    .. math::

        \vec{p} = -\nu \frac{S}{c} \left ( \frac{C_{r}A}{m} \right )\frac{\vec{r}}{r}

    Parameters
    ----------
    t0 : float
        Current time (s)
    state : numpy.ndarray
        Six component state vector [x, y, z, vx, vy, vz] (km, km/s).
    k : float
        gravitational constant, (km^3/s^2)
    R : float
        radius of the attractor
    C_R: float
        dimensionless radiation pressure coefficient, 1 < C_R < 2 ()
    A_over_m: float
        effective spacecraft area/mass of the spacecraft (km^2/kg)
    Wdivc_s : float
        total star emitted power divided by the speed of light (W * s / km)
    star: a callable object returning the position of star in attractor frame
        star position

    Note
    ----
    This function provides the acceleration due to star light pressure. We follow
    Howard Curtis, section 12.9

    """

    r_star = star(t0)
    r_sat = state[:3]
    P_s = Wdivc_s / (norm(r_star) ** 2)

    nu = float(shadow_function(r_sat, r_star, R))
    return -nu * P_s * (C_R * A_over_m) * r_star / norm(r_star)
