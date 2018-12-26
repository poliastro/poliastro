import numpy as np
from astropy import units as u

from poliastro.util import norm


@u.quantity_input(
    v_spacecraft=u.km / u.s,
    v_body=u.km / u.s,
    k=u.km ** 3 / u.s ** 2,
    r_p=u.km,
    theta=u.deg,
)
def compute_flyby(v_spacecraft, v_body, k, r_p, theta=0 * u.deg):
    """Computes outbound velocity after a flyby.

    Parameters
    ----------
    v_spacecraft : ~astropy.units.Quantity
        Velocity of the spacecraft, relative to the attractor of the body.
    v_body : ~astropy.units.Quantity
        Velocity of the body, relative to its attractor.
    k : ~astropy.units.Quantity
        Standard gravitational parameter of the body.
    r_p : ~astropy.units.Quantity
        Radius of periapsis, measured from the center of the body.
    theta : ~astropy.units.Quantity, optional
        Aim angle of the B vector, default to 0.

    Returns
    -------
    v_spacecraft_out : ~astropy.units.Quantity
        Outbound velocity of the spacecraft.
    delta : ~astropy.units.Quantity
        Turn angle.

    """
    v_inf_1 = v_spacecraft - v_body  # Hyperbolic excess velocity
    v_inf = norm(v_inf_1)

    ecc = 1 + r_p * v_inf ** 2 / k  # Eccentricity of the entry hyperbola
    delta = 2 * np.arcsin(1 / ecc)  # Turn angle

    b = k / v_inf ** 2 * np.sqrt(ecc ** 2 - 1)  # Magnitude of the B vector

    # Now we compute the unit vectors in which to return the outbound hyperbolic excess velocity:
    # * S goes along the hyperbolic excess velocity and is perpendicular to the B-Plane,
    # * T goes along the B-Plane and is parallel to _some_ fundamental plane - in this case, the plane in which
    #   the velocities are computed
    # * R completes the orthonormal set
    S_vec = v_inf_1 / v_inf
    c_vec = np.array([0, 0, 1]) * u.one
    T_vec = np.cross(S_vec, c_vec) * u.one
    T_vec = T_vec / norm(T_vec)
    R_vec = np.cross(S_vec, T_vec) * u.one

    # This vector defines the B-Plane
    B_vec = b * (np.cos(theta) * T_vec + np.sin(theta) * R_vec)

    # We have to rotate the inbound hyperbolic excess velocity
    # an angle delta (turn angle) around a vector that is orthogonal to
    # the B-Plane and trajectory plane
    rot_v = np.cross(B_vec / b, S_vec) * u.one

    # And now we rotate the outbound hyperbolic excess velocity
    # u_vec = v_inf_1 / norm(v_inf) = S_vec
    v_vec = np.cross(rot_v, v_inf_1) * u.one
    v_vec = v_vec / norm(v_vec)

    v_inf_2 = v_inf * (np.cos(delta) * S_vec + np.sin(delta) * v_vec)

    v_spacecraft_out = v_inf_2 + v_body

    return v_spacecraft_out, delta
