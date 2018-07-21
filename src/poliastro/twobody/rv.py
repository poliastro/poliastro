"""Functions to define orbits from position and velocity vectors.

"""
from astropy import units as u

from poliastro.twobody import classical
from poliastro.core.elements import rv2coe

from ._base import BaseState


class RVState(BaseState):
    """State defined by its position and velocity vectors.

    """
    def __init__(self, attractor, r, v):
        super().__init__(attractor)
        self._r = r
        self._v = v

    @property
    def r(self):
        """Position vector. """
        return self._r

    @property
    def v(self):
        """Velocity vector. """
        return self._v

    def to_vectors(self):
        """Converts to position and velocity vector representation.

        """
        return self

    def to_classical(self):
        """Converts to classical orbital elements representation.

        """
        (p, ecc, inc, raan, argp, nu
         ) = rv2coe(self.attractor.k.to(u.km ** 3 / u.s ** 2).value,
                    self.r.to(u.km).value,
                    self.v.to(u.km / u.s).value)

        return classical.ClassicalState(
            self.attractor,
            p * u.km,
            ecc * u.one,
            inc * u.rad,
            raan * u.rad,
            argp * u.rad,
            nu * u.rad)
