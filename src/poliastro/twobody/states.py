from astropy import units as u

from poliastro.core.elements import coe2mee, coe2rv, mee2coe, mee2rv, rv2coe
from poliastro.twobody.elements import mean_motion, period


class BaseState:
    """Base State class, meant to be subclassed."""

    def __init__(self, attractor, plane):
        """Constructor.

        Parameters
        ----------
        attractor : Body
            Main attractor.
        plane : ~poliastro.frames.enums.Planes
            Reference plane for the elements.

        """
        self._attractor = attractor
        self._plane = plane

    @property
    def plane(self):
        """Fundamental plane of the frame."""
        return self._plane

    @property
    def attractor(self):
        """Main attractor."""
        return self._attractor

    @property
    def n(self):
        """Mean motion."""
        return mean_motion(self.attractor.k, self.to_classical().a)

    @property
    def period(self):
        """Period of the orbit."""
        return period(self.attractor.k, self.to_classical().a)

    def to_vectors(self):
        """Converts to position and velocity vector representation.

        Returns
        -------
        RVState

        """
        raise NotImplementedError

    def to_classical(self):
        """Converts to classical orbital elements representation.

        Returns
        -------
        ClassicalState

        """
        raise NotImplementedError

    def to_equinoctial(self):
        """Converts to modified equinoctial elements representation.

        Returns
        -------
        ModifiedEquinoctialState

        """
        raise NotImplementedError


class ClassicalState(BaseState):
    """State defined by its classical orbital elements."""

    def __init__(self, attractor, p, ecc, inc, raan, argp, nu, plane):
        """Constructor.

        Parameters
        ----------
        attractor : Body
            Main attractor.
        p : ~astropy.units.Quantity
            Semilatus rectum.
        ecc : ~astropy.units.Quantity
            Eccentricity.
        inc : ~astropy.units.Quantity
            Inclination.
        raan : ~astropy.units.Quantity
            Right ascension of the ascending node.
        argp : ~astropy.units.Quantity
            Argument of the perigee.
        nu : ~astropy.units.Quantity
            True anomaly.
        plane : ~poliastro.frames.enums.Planes
            Reference plane for the elements.

        """
        super().__init__(attractor, plane)
        self._p = p
        self._ecc = ecc
        self._inc = inc
        self._raan = raan
        self._argp = argp
        self._nu = nu

    @property
    def p(self):
        """Semilatus rectum."""
        return self._p

    @property
    def a(self):
        """Semimajor axis."""
        return self.p / (1 - self.ecc ** 2)

    @property
    def ecc(self):
        """Eccentricity."""
        return self._ecc

    @property
    def inc(self):
        """Inclination."""
        return self._inc

    @property
    def raan(self):
        """Right ascension of the ascending node."""
        return self._raan

    @property
    def argp(self):
        """Argument of the perigee."""
        return self._argp

    @property
    def nu(self):
        """True anomaly."""
        return self._nu

    def to_vectors(self):
        """Converts to position and velocity vector representation."""
        r, v = coe2rv(
            self.attractor.k.to_value(u.km ** 3 / u.s ** 2),
            self.p.to_value(u.km),
            self.ecc.value,
            self.inc.to_value(u.rad),
            self.raan.to_value(u.rad),
            self.argp.to_value(u.rad),
            self.nu.to_value(u.rad),
        )

        return RVState(self.attractor, r * u.km, v * u.km / u.s, self.plane)

    def to_classical(self):
        """Converts to classical orbital elements representation."""
        return self

    def to_equinoctial(self):
        """Converts to modified equinoctial elements representation."""
        p, f, g, h, k, L = coe2mee(
            self.p.to_value(u.km),
            self.ecc.value,
            self.inc.to_value(u.rad),
            self.raan.to_value(u.rad),
            self.argp.to_value(u.rad),
            self.nu.to_value(u.rad),
        )

        return ModifiedEquinoctialState(
            self.attractor,
            p * u.km,
            f * u.rad,
            g * u.rad,
            h * u.rad,
            k * u.rad,
            L * u.rad,
            self.plane,
        )


class RVState(BaseState):
    """State defined by its position and velocity vectors."""

    def __init__(self, attractor, r, v, plane):
        """Constructor.

        Parameters
        ----------
        attractor : Body
            Main attractor.
        r : ~astropy.units.Quantity
            Position vector wrt attractor center.
        v : ~astropy.units.Quantity
            Velocity vector.
        plane : ~poliastro.frames.enums.Planes
            Reference plane for the elements.

        """
        super().__init__(attractor, plane)
        self._r = r
        self._v = v

    @property
    def r(self):
        """Position vector."""
        return self._r

    @property
    def v(self):
        """Velocity vector."""
        return self._v

    def to_vectors(self):
        """Converts to position and velocity vector representation."""
        return self

    def to_classical(self):
        """Converts to classical orbital elements representation."""
        (p, ecc, inc, raan, argp, nu) = rv2coe(
            self.attractor.k.to_value(u.km ** 3 / u.s ** 2),
            self.r.to_value(u.km),
            self.v.to_value(u.km / u.s),
        )

        return ClassicalState(
            self.attractor,
            p * u.km,
            ecc * u.one,
            inc * u.rad,
            raan * u.rad,
            argp * u.rad,
            nu * u.rad,
            self.plane,
        )


class ModifiedEquinoctialState(BaseState):
    """State defined by modified equinoctial elements representation."""

    def __init__(self, attractor, p, f, g, h, k, L, plane):
        super().__init__(attractor, plane)
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
            self.p.to_value(u.km),
            self.f.to_value(u.rad),
            self.g.to_value(u.rad),
            self.h.to_value(u.rad),
            self.k.to_value(u.rad),
            self.L.to_value(u.rad),
        )

        return ClassicalState(
            self.attractor,
            p * u.km,
            ecc * u.one,
            inc * u.rad,
            raan * u.rad,
            argp * u.rad,
            nu * u.rad,
            self.plane,
        )

    def to_vectors(self):
        r, v = mee2rv(
            self.p.to(u.km).value,
            self.f.to(u.rad).value,
            self.g.to(u.rad).value,
            self.h.to(u.rad).value,
            self.k.to(u.rad).value,
            self.L.to(u.rad).value,
        )
        return RVState(self.attractor, r * u.km, v * u.km / u.s, self.plane)
