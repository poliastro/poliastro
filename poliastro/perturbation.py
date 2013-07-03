"""Perturbation analysis.

"""

import numpy as np
from numpy import sin, cos, tan, sqrt

from poliastro import angles
from poliastro.logging import logger


def gauss(x, t, k, fp, kw_fp):
    """Gauss' form of the Lagrange plantary equations.

    Return the derivative of the elements, given the vector of elements,
    the time and the gravitational parameter of the attractor.

    Parameters
    ----------
    x : array_like
        Vector of classical orbital elements: a, ecc, inc, omega, argp, M.
    t : float
        Seconds elapsed from EPOCH.
    k : float
        Gravitational parameter.
    fp : function
        Perturbing force per unit mass function.
    kw_fp : dict
        Extra keyword arguments of the perturbing function.

    Notes
    -----
    These are the VOP equations in Gauss form, as derived in [1]_ and
    conveniently modified. They are expressed in the RSW system.

    References
    ----------
    .. [1] Bate, "Fundamentals of Astrodynamics", Dover, pp. 396-406, 1971.

    """
    logger.debug(t)
    a, ecc, inc, omega, argp, M = x
    p = a * (1 - ecc ** 2)
    n = sqrt(k / a ** 3)
    _, nu = angles.M2nu(ecc, M)
    u = argp + nu
    coe = p, ecc, inc, omega, argp, nu
    f_R, f_S, f_W = fp(t, k, *coe, **kw_fp)
    h = sqrt(k * p)
    r = p / (1 + ecc * cos(nu))
    dx = np.array([
        2 * a ** 2 / h * (f_R * ecc * sin(nu) + p / r * f_S),
        1 / h * (p * sin(nu) * f_R + ((p + r) * cos(nu) + r * ecc) * f_S),
        1 / h * cos(u) * r * f_W,
        1 / h * sin(u) / sin(inc) * r * f_W,
        (1 / (ecc * h) * (-p * cos(nu) * f_R + (p + r) * sin(nu) * f_S) -
         1 / h * sin(u) / tan(inc) * r * f_W),
        (n + sqrt(1 - ecc ** 2) / (ecc * h) *
         ((p * cos(nu) - 2 * r * ecc) * f_R - (p + r) * sin(nu) * f_S))
    ])
    return dx
