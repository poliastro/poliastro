import numpy as np
from numba import njit as jit

from poliastro._math.linalg import norm
from poliastro.core.elements import coe_rotation_matrix, rv2coe
from poliastro.core.util import planetocentric_to_AltAz


@jit
def eclipse_function(k, u_, r_sec, R_sec, R_primary, umbra=True):
    """Calculates a continuous shadow function.

    Parameters
    ----------
    k : float
        Standard gravitational parameter (km^3 / s^2).
    u_ : numpy.ndarray
        Satellite position and velocity vector with respect to the primary body.
    r_sec : numpy.ndarray
        Position vector of the secondary body with respect to the primary body.
    R_sec : float
        Equatorial radius of the secondary body.
    R_primary : float
        Equatorial radius of the primary body.
    umbra : bool
        Whether to calculate the shadow function for umbra or penumbra, defaults to True
        i.e. calculates for umbra.

    Notes
    -----
    The shadow function is taken from Escobal, P. (1985). Methods of orbit determination.
    The current implementation assumes circular bodies and doesn't account for flattening.

    """
    # Plus or minus condition
    pm = 1 if umbra else -1
    p, ecc, inc, raan, argp, nu = rv2coe(k, u_[:3], u_[3:])

    PQW = coe_rotation_matrix(inc, raan, argp)
    # Make arrays contiguous for faster dot product with numba.
    P_, Q_ = np.ascontiguousarray(PQW[:, 0]), np.ascontiguousarray(PQW[:, 1])

    r_sec_norm = norm(r_sec)
    beta = (P_ @ r_sec) / r_sec_norm
    zeta = (Q_ @ r_sec) / r_sec_norm

    sin_delta_shadow = np.sin((R_sec - pm * R_primary) / r_sec_norm)

    cos_psi = beta * np.cos(nu) + zeta * np.sin(nu)
    shadow_function = (
        ((R_primary**2) * (1 + ecc * np.cos(nu)) ** 2)
        + (p**2) * (cos_psi**2)
        - p**2
        + pm
        * (2 * p * R_primary * cos_psi)
        * (1 + ecc * np.cos(nu))
        * sin_delta_shadow
    )

    return shadow_function


@jit
def line_of_sight(r1, r2, R):
    """Calculates the line of sight condition between two position vectors, r1 and r2.

    Parameters
    ----------
    r1 : numpy.ndarray
        The position vector of the first object with respect to a central attractor.
    r2 : numpy.ndarray
        The position vector of the second object with respect to a central attractor.
    R : float
        The radius of the central attractor.

    Returns
    -------
    delta_theta: float
        Greater than or equal to zero, if there exists a LOS between two objects
        located by r1 and r2, else negative.

    """
    r1_norm = norm(r1)
    r2_norm = norm(r2)

    theta = np.arccos((r1 @ r2) / r1_norm / r2_norm)
    theta_1 = np.arccos(R / r1_norm)
    theta_2 = np.arccos(R / r2_norm)

    return (theta_1 + theta_2) - theta


@jit
def elevation_function(k, u_, phi, theta, R, R_p, H):
    """Calculates the elevation angle of an object in orbit with respect to
    a location on attractor.

    Parameters
    ----------
    k: float
        Standard gravitational parameter.
    u_: numpy.ndarray
        Satellite position and velocity vector with respect to the central attractor.
    phi: float
        Geodetic Latitude of the station.
    theta: float
        Local sidereal time at a particular instant.
    R: float
        Equatorial radius of the central attractor.
    R_p: float
        Polar radius of the central attractor.
    H: float
        Elevation, above the ellipsoidal surface.
    """
    ecc = np.sqrt(1 - (R_p / R) ** 2)
    denom = np.sqrt(1 - ecc**2 * np.sin(phi) ** 2)
    g1 = H + (R / denom)
    g2 = H + (1 - ecc**2) * R / denom
    # Coordinates of location on attractor.
    coords = np.array(
        [
            g1 * np.cos(phi) * np.cos(theta),
            g1 * np.cos(phi) * np.sin(theta),
            g2 * np.sin(phi),
        ]
    )

    # Position of satellite with respect to a point on attractor.
    rho = np.subtract(u_[:3], coords)

    rot_matrix = planetocentric_to_AltAz(theta, phi)

    new_rho = rot_matrix @ rho
    new_rho = new_rho / np.linalg.norm(new_rho)
    el = np.arcsin(new_rho[-1])

    return el
