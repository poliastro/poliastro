import numpy as np
from astropy import units as u

from ..util import alinspace
from .angles import E_to_nu, nu_to_E


@u.quantity_input(min_nu=u.rad, ecc=u.one, max_nu=u.rad)
def sample_closed(min_nu, ecc, max_nu=None, num_values=100):
    """Sample a closed orbit

    If ``max_nu`` is given, the sampling interval will go
    from the minimum to the maximum true anomaly in the direction of the orbit.
    If not given, it will do a full revolution starting in the minimum true anomaly.

    Notes
    -----
    First sample the eccentric anomaly uniformly,
    then transform into true anomaly
    to minimize error in the apocenter,
    see https://apps.dtic.mil/dtic/tr/fulltext/u2/a605040.pdf

    """
    # Because how nu_to_E works, we don't need to wrap the angle here!
    # It will do the right thing
    min_E = nu_to_E(min_nu, ecc)

    # This linspace will always increase positively,
    # even though it might contain out of range values
    E_values = alinspace(
        min_E, nu_to_E(max_nu, ecc) if max_nu is not None else None, num=num_values
    )

    # Because how E_to_nu works, we don't need to wrap the angles here!
    # It will do the right thing
    nu_values = E_to_nu(E_values, ecc)

    # We wrap the angles on return
    return (nu_values + np.pi * u.rad) % (2 * np.pi * u.rad) - np.pi * u.rad
