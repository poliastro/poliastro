"""Functions to define orbits from modified equinoctial orbital elements.

"""
from astropy import units as u

from poliastro.core.elements import mee2coe
from poliastro.twobody import classical

from ._base import BaseState


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
        p, ecc, inc, raan, argp, nu = mee2coe(
            self.p.to(u.km).value,
            self.f.to(u.rad).value,
            self.g.to(u.rad).value,
            self.h.to(u.rad).value,
            self.k.to(u.rad).value,
            self.L.to(u.rad).value,
        )

        return classical.ClassicalState(
            self.attractor,
            p * u.km,
            ecc * u.one,
            inc * u.rad,
            raan * u.rad,
            argp * u.rad,
            nu * u.rad,
        )
