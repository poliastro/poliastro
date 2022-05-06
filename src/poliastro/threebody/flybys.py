from astropy import units as u

from poliastro.core.flybys import compute_flyby as compute_flyby_fast


@u.quantity_input(
    v_spacecraft=u.km / u.s,
    v_body=u.km / u.s,
    k=u.km**3 / u.s**2,
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
    v_spacecraft = v_spacecraft.to_value(u.km / u.s)
    v_body = v_body.to_value(u.km / u.s)
    k = k.to_value(u.km**3 / u.s**2)
    r_p = r_p.to_value(u.km)
    theta = theta.to_value(u.rad)

    v_spacecraft_out, delta = compute_flyby_fast(
        v_spacecraft, v_body, k, r_p, theta
    )

    return v_spacecraft_out * u.km / u.s, delta * u.rad
