"""Functions to define orbits from classical orbital elements.

"""
from astropy import units as u

from poliastro.core.elements import coe2mee, coe2rv
from poliastro.twobody import equinoctial, rv

from ._base import BaseState


class ClassicalState(BaseState):
    """State defined by its classical orbital elements.

    """

    def __init__(self, attractor, p, ecc, inc, raan, argp, nu):
        super().__init__(attractor)
        self._p = p
        self._ecc = ecc
        self._inc = inc
        self._raan = raan
        self._argp = argp
        self._nu = nu

    @property
    def p(self):
        """Semilatus rectum. """
        return self._p

    @property
    def ecc(self):
        """Eccentricity. """
        return self._ecc

    @property
    def inc(self):
        """Inclination. """
        return self._inc

    @property
    def raan(self):
        """Right ascension of the ascending node. """
        return self._raan

    @property
    def argp(self):
        """Argument of the perigee. """
        return self._argp

    @property
    def nu(self):
        """True anomaly. """
        return self._nu

    def to_vectors(self):
        """Converts to position and velocity vector representation.

        """
        r, v = coe2rv(
            self.attractor.k.to(u.km ** 3 / u.s ** 2).value,
            self.p.to(u.km).value,
            self.ecc.value,
            self.inc.to(u.rad).value,
            self.raan.to(u.rad).value,
            self.argp.to(u.rad).value,
            self.nu.to(u.rad).value,
        )

        return rv.RVState(self.attractor, r * u.km, v * u.km / u.s)

    def to_classical(self):
        """Converts to classical orbital elements representation.

        """
        return self

    def to_equinoctial(self):
        """Converts to modified equinoctial elements representation.

        """
        p, f, g, h, k, L = coe2mee(
            self.p.to(u.km).value,
            self.ecc.value,
            self.inc.to(u.rad).value,
            self.raan.to(u.rad).value,
            self.argp.to(u.rad).value,
            self.nu.to(u.rad).value,
        )

        return equinoctial.ModifiedEquinoctialState(
            self.attractor,
            p * u.km,
            f * u.rad,
            g * u.rad,
            h * u.rad,
            k * u.rad,
            L * u.rad,
        )
