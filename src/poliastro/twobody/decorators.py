"""Decorators.

"""

from functools import wraps

from astropy import units as u

from poliastro.bodies import Body
from poliastro.twobody.rv import RVState

u.kms = u.km / u.s
u.km3s2 = u.km ** 3 / u.s ** 2


def state_from_vector(func):

    @wraps(func)
    def wrapper(t, u_, k, *args, **kwargs):
        r, v = u_[:3], u_[3:]
        ss = RVState(Body(None, k * u.km3s2, "_Dummy"), r * u.km, v * u.kms)

        return func(t, ss, *args, **kwargs)

    return wrapper
