"""Functions to define orbits from modified equinoctial orbital elements.

"""
import numpy as np
from astropy import units as u

from poliastro.twobody import classical

from ._base import BaseState


def mee2coe(p, f, g, h, k, L):
    """Converts from modified equinoctial orbital elements to classical
    orbital elements.

    The definition of the modified equinoctial orbital elements is taken from
    [Walker, 1985].

    Note
    -----
    The conversion is always safe because arctan2 works also for 0, 0
    arguments.

    """
    ecc = np.sqrt(f**2 + g**2)
    inc = 2 * np.arctan(np.sqrt(h**2 + k**2))
    lonper = np.arctan2(g, f)
    raan = np.arctan2(k, h) % (2 * np.pi)
    argp = (lonper - raan) % (2 * np.pi)
    nu = (L - lonper) % (2 * np.pi)
    return p, ecc, inc, raan, argp, nu


class ModifiedEquinoctialState(BaseState):
    def __init__(self, attractor, p, f, g, h, k, L):
        super().__init__(attractor)
        self._p = p
        self._f = f
        self._g = g
        self._h = h
        self._k = k
        self._L = L

    @property
    def p(self):
        return self._p

    @property
    def f(self):
        return self._f

    @property
    def g(self):
        return self._g

    @property
    def h(self):
        return self._h

    @property
    def k(self):
        return self._k

    @property
    def L(self):
        return self._L

    def to_classical(self):
        p, ecc, inc, raan, argp, nu = mee2coe(self.p.to(u.km).value,
                                              self.f.to(u.rad).value,
                                              self.g.to(u.rad).value,
                                              self.h.to(u.rad).value,
                                              self.k.to(u.rad).value,
                                              self.L.to(u.rad).value)

        return classical.ClassicalState(
            self.attractor,
            p * u.km,
            ecc * u.one,
            inc * u.rad,
            raan * u.rad,
            argp * u.rad,
            nu * u.rad)
