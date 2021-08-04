import numpy as np
from numba import njit as jit
from numpy.linalg import norm

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
def satellite_visibility(u_, N, phi, H, R_primary):
    p, ecc, inc, raan, argp, nu = rv2coe(k, u_[:3], u_[3:])
    a = p / (1 - ecc ** 2)

    PQW = coe_rotation_matrix(inc, raan, argp)
    P_, Q_ = np.ascontiguousarray(PQW[:, 0]), np.ascontiguousarray(PQW[:, 1])

    denom = np.sqrt(1 - (2 * f - f ** 2) * np.sin(phi) ** 2))
    g1 = H + (R_primary / denom)
    g2 = H + (1 - f ** 2) * R_primary / denom
    g = g1 * np.cos(phi) ** 2 + g2 * np.sin(phi) ** 2
    visibility_function = a * (np.cos(E) - ecc) * np.dot(P_, N) + (a * np.sqrt(1 - ecc ** 2) * np.sin(E)) * np.dot(Q_, N) - g

    return visibility_function
