""" Low level computations for flybys """

import numpy as np
from numba import njit as jit
from numpy import cross

from poliastro._math.linalg import norm


@jit
def compute_flyby(v_spacecraft, v_body, k, r_p, theta):
    """Computes outbound velocity after a flyby and the turn angle

    Parameters
    ----------
    k : float
        Standard Gravitational parameter.
    r_p : float
        Radius of periapsis, measured from the center of the body.
    v_spacecraft : float
        Velocity of the spacecraft, relative to the attractor of the body.
    v_body : float
        Velocity of the body, relative to its attractor.
    theta : float
        Aim angle of the B vector.

    Returns
    -------
    v_spacecraft_out : float
        Outbound velocity of the spacecraft.
    delta : float
        Turn angle.

    """
    v_inf_1 = v_spacecraft - v_body  # Hyperbolic excess velocity
    v_inf = norm(v_inf_1)

    ecc = 1 + r_p * v_inf**2 / k  # Eccentricity of the entry hyperbola
    delta = 2 * np.arcsin(1 / ecc)  # Turn angle

    b = k / v_inf**2 * np.sqrt(ecc**2 - 1)  # Magnitude of the B vector

    # Now we compute the unit vectors in which to return the outbound hyperbolic excess velocity:
    # * S goes along the hyperbolic excess velocity and is perpendicular to the B-Plane,
    # * T goes along the B-Plane and is parallel to _some_ fundamental plane - in this case, the plane in which
    #   the velocities are computed
    # * R completes the orthonormal set
    S_vec = v_inf_1 / v_inf
    c_vec = np.array([0, 0, 1])
    T_vec = cross(S_vec, c_vec)
    T_vec = T_vec / norm(T_vec)
    R_vec = cross(S_vec, T_vec)

    # This vector defines the B-Plane
    B_vec = b * (np.cos(theta) * T_vec + np.sin(theta) * R_vec)

    # We have to rotate the inbound hyperbolic excess velocity
    # an angle delta (turn angle) around a vector that is orthogonal to
    # the B-Plane and trajectory plane
    rot_v = cross(B_vec / b, S_vec)

    # And now we rotate the outbound hyperbolic excess velocity
    # u_vec = v_inf_1 / norm(v_inf) = S_vec
    v_vec = cross(rot_v, v_inf_1)
    v_vec = v_vec / norm(v_vec)

    v_inf_2 = v_inf * (np.cos(delta) * S_vec + np.sin(delta) * v_vec)

    v_spacecraft_out = v_inf_2 + v_body

    return v_spacecraft_out, delta
