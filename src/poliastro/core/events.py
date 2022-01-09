import numpy as np
from numba import njit as jit
from numpy.linalg import norm

from poliastro.core.elements import coe_rotation_matrix, rv2coe


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
    r1_norm = np.linalg.norm(r1)
    r2_norm = np.linalg.norm(r2)

    theta = np.arccos((r1 @ r2) / r1_norm / r2_norm)
    theta_1 = np.arccos(R / r1_norm)
    theta_2 = np.arccos(R / r2_norm)

    return (theta_1 + theta_2) - theta
