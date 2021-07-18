import numpy as np
from numba import njit as jit
from numpy.linalg import norm

from poliastro.core.elements import rv2coe


@jit
def continuous_shadow_function_for_penumbra(k, u_, r_sec, R_sec, R_primary):
    p, ecc, inc, raan, argp, nu = rv2coe(k, u_[:3], u_[3:])

    Px = np.cos(raan) * np.cos(argp) - np.sin(raan) * np.sin(argp) * np.cos(inc)
    Py = np.sin(raan) * np.cos(argp) + np.cos(raan) * np.sin(argp) * np.cos(inc)
    Pz = np.sin(argp) * np.sin(inc)
    P_ = np.array([Px, Py, Pz])
    Qx = -np.cos(raan) * np.sin(argp) - np.sin(raan) * np.cos(argp) * np.cos(inc)
    Qy = -np.sin(raan) * np.sin(argp) + np.cos(raan) * np.cos(argp) * np.cos(inc)
    Qz = np.cos(argp) * np.sin(inc)
    Q_ = np.array([Qx, Qy, Qz])

    r_sec_norm = norm(r_sec)
    beta = np.dot(P_, r_sec) / r_sec_norm
    zeta = np.dot(Q_, r_sec) / r_sec_norm

    delta_p = np.arcsin((R_sec + R_primary) / r_sec_norm)

    cos_psi = beta * np.cos(nu) + zeta * np.sin(nu)
    shadow_function = (
        ((R_primary ** 2) * (1 + ecc * np.cos(nu)) ** 2)
        + (p ** 2) * (cos_psi ** 2)
        - p ** 2
        - (2 * p * R_primary * cos_psi) * (1 + ecc * np.cos(nu)) * np.sin(delta_p)
    )

    return shadow_function


@jit
def continuous_shadow_function_for_umbra(k, u_, r_sec, R_sec, R_primary):
    p, ecc, inc, raan, argp, nu = rv2coe(k, u_[:3], u_[3:])

    Px = np.cos(raan) * np.cos(argp) - np.sin(raan) * np.sin(argp) * np.cos(inc)
    Py = np.sin(raan) * np.cos(argp) + np.cos(raan) * np.sin(argp) * np.cos(inc)
    Pz = np.sin(argp) * np.sin(inc)
    P_ = np.array([Px, Py, Pz])
    Qx = -np.cos(raan) * np.sin(argp) - np.sin(raan) * np.cos(argp) * np.cos(inc)
    Qy = -np.sin(raan) * np.sin(argp) + np.cos(raan) * np.cos(argp) * np.cos(inc)
    Qz = np.cos(argp) * np.sin(inc)
    Q_ = np.array([Qx, Qy, Qz])

    r_sec_norm = norm(r_sec)
    beta = np.dot(P_, r_sec) / r_sec_norm
    zeta = np.dot(Q_, r_sec) / r_sec_norm

    delta_u = np.arcsin((R_sec - R_primary) / r_sec_norm)

    cos_psi = beta * np.cos(nu) + zeta * np.sin(nu)
    shadow_function = (
        ((R_primary ** 2) * (1 + ecc * np.cos(nu)) ** 2)
        + (p ** 2) * (cos_psi ** 2)
        - p ** 2
        + (2 * p * R_primary * cos_psi) * (1 + ecc * np.cos(nu)) * np.sin(delta_u)
    )

    return shadow_function
