"""Decorators.

"""
from functools import wraps

from astropy import units as u

from poliastro.bodies import Body
from poliastro.frames import Planes
from poliastro.twobody.states import RVState

u.kms = u.km / u.s
u.km3s2 = u.km**3 / u.s**2


def state_from_vector(func):
    """Changes signature to receive Orbit instead of state array.

    Examples
    --------
    >>> from poliastro.twobody.decorators import state_from_vector
    >>> @state_from_vector
    ... def func(_, ss):
    ...     return ss.r, ss.v
    ...
    >>> func(0.0, [1, 2, 3, -1, -2, -3], 1.0)
    (<Quantity [1., 2., 3.] km>, <Quantity [-1., -2., -3.] km / s>)

    Notes
    -----
    Functions decorated with this will have poor performance.

    """

    @wraps(func)
    def wrapper(t, u_, k, *args, **kwargs):
        r, v = u_[:3], u_[3:]
        ss = RVState(
            Body(None, k * u.km3s2, "_Dummy"),
            (
                r * u.km,
                v * u.kms,
            ),
            Planes.EARTH_EQUATOR,
        )

        return func(t, ss, *args, **kwargs)

    return wrapper
