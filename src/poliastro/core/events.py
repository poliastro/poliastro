import numpy as np
from numba import njit as jit
from numpy.linalg import norm

from poliastro.core.angles import nu_to_E
from poliastro.core.elements import coe_rotation_matrix, rv2coe


@jit
def eclipse_function(k, u_, r_sec, R_sec, R_primary, umbra=True):
    """Calculates a continuous shadow function.

    Parameters
    ----------
    k: float
        Standard gravitational parameter (km^3 / s^2).
    u_: ~np.array
        Satellite position and velocity vector with respect to the primary body.
    r_sec: ~np.array
        Position vector of the secondary body with respect to the primary body.
    R_sec: float
        Equatorial radius of the secondary body.
    R_primary: float
        Equatorial radius of the primary body.
    umbra: bool
        Whether to calculate the shadow function for umbra or penumbra, defaults to True
        i.e. calculates for umbra.

    Note
    ----
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
    beta = np.dot(P_, r_sec) / r_sec_norm
    zeta = np.dot(Q_, r_sec) / r_sec_norm

    sin_delta_shadow = np.sin((R_sec - pm * R_primary) / r_sec_norm)

    cos_psi = beta * np.cos(nu) + zeta * np.sin(nu)
    shadow_function = (
        ((R_primary ** 2) * (1 + ecc * np.cos(nu)) ** 2)
        + (p ** 2) * (cos_psi ** 2)
        - p ** 2
        + pm * (2 * p * R_primary * cos_psi) * (1 + ecc * np.cos(nu)) * sin_delta_shadow
    )

    return shadow_function


@jit
def line_of_sight(r1, r2, R):
    """Calculates the line of sight condition between two position vectors, r1 and r2.

    Parameters
    ----------
    r1: ~np.array
        The position vector of the first object with respect to a central attractor.
    r2: ~np.array
        The position vector of the second object with respect to a central attractor.
    R: float
        The radius of the central attractor.

    Returns
    -------
    delta_theta: float
        Greater than or equal to zero, if there exists a LOS between two objects
        located by r1 and r2, else negative.

    """
    r1_norm = np.linalg.norm(r1)
    r2_norm = np.linalg.norm(r2)

    theta = np.arccos(np.dot(r1, r2) / r1_norm / r2_norm)
    theta_1 = np.arccos(R / r1_norm)
    theta_2 = np.arccos(R / r2_norm)

    return (theta_1 + theta_2) - theta


@jit
def satellite_view(k, u_, r_sec):
    """Calculates a continuous visibility function to check whether a body with a
    given position vector, u_[:3], can view the attractor's surface illuminated by the secondary
    body, with position vector `r_sec`. All the position vectors are with respect to the attractor.

    Parameters
    ----------
    k: float
        Standard gravitational parameter (km^3 / s^2).
    u_: ~numpy.array
        Satellite position and velocity vector with respect to the primary body
        i.e the body around which the satellite revolves.
    r_sec: ~numpy.array
        Position vector of the secondary body with respect to the primary body.

    Note
    ----
    The implementation supports only elliptical orbits.

    """
    p, ecc, inc, raan, argp, nu = rv2coe(k, u_[:3], u_[3:])
    E = nu_to_E(nu, ecc)

    r_sec_norm = norm(r_sec)
    r_norm = norm(u_[:3])

    PQW = coe_rotation_matrix(inc, raan, argp)
    P_, Q_ = np.ascontiguousarray(PQW[:, 0]), np.ascontiguousarray(PQW[:, 1])

    cos_psi = np.dot(u_[:3], r_sec) / r_sec_norm / r_norm
    view_function = (
        (np.cos(E) - ecc) * (np.dot(P_, r_sec))
        + (np.sqrt(1 - ecc ** 2) * np.sin(E)) * np.dot(Q_, r_sec)
        - (1 - ecc * np.cos(E)) * r_sec_norm * cos_psi
    )

    return view_function
