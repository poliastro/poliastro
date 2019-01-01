import numpy as np
from astropy import units as u

from poliastro.core.elements import coe2mee, coe2rv, mee2coe, rv2coe
from poliastro.util import norm


class BaseState(object):
    """Base State class, meant to be subclassed.

    """

    def __init__(self, attractor):
        """Constructor.

        Parameters
        ----------
        attractor : Body
            Main attractor.

        """
        self._attractor = attractor

    @property
    def attractor(self):
        """Main attractor. """
        return self._attractor

    @property
    def r(self):
        """Position vector. """
        return self.to_vectors().r

    @property
    def v(self):
        """Velocity vector. """
        return self.to_vectors().v

    @property
    def a(self):
        """Semilatus rectum. """
        return self.p / (1 - self.ecc ** 2)

    @property
    def p(self):
        """Semimajor axis. """
        return self.to_classical().p

    @property
    def r_p(self):
        """Radius of pericenter. """
        return self.a * (1 - self.ecc)

    @property
    def r_a(self):
        """Radius of apocenter. """
        return self.a * (1 + self.ecc)

    @property
    def ecc(self):
        """Eccentricity. """
        return self.to_classical().ecc

    @property
    def inc(self):
        """Inclination. """
        return self.to_classical().inc

    @property
    def raan(self):
        """Right ascension of the ascending node. """
        return self.to_classical().raan

    @property
    def argp(self):
        """Argument of the perigee. """
        return self.to_classical().argp

    @property
    def nu(self):
        """True anomaly. """
        return self.to_classical().nu

    @property
    def f(self):
        """Second modified equinoctial element. """
        return self.to_equinoctial().f

    @property
    def g(self):
        """Third modified equinoctial element. """
        return self.to_equinoctial().g

    @property
    def h(self):
        """Fourth modified equinoctial element. """
        return self.to_equinoctial().h

    @property
    def k(self):
        """Fifth modified equinoctial element. """
        return self.to_equinoctial().k

    @property
    def L(self):
        """True longitude. """
        return self.raan + self.argp + self.nu

    @property
    def period(self):
        """Period of the orbit. """
        return 2 * np.pi * u.rad / self.n

    @property
    def n(self):
        """Mean motion. """
        return (np.sqrt(self.attractor.k / abs(self.a ** 3)) * u.rad).decompose()

    @property
    def energy(self):
        """Specific energy. """
        return self.v.dot(self.v) / 2 - self.attractor.k / np.sqrt(self.r.dot(self.r))

    @property
    def e_vec(self):
        """Eccentricity vector. """
        r, v = self.rv()
        k = self.attractor.k
        e_vec = ((v.dot(v) - k / (norm(r))) * r - r.dot(v) * v) / k
        return e_vec.decompose()

    @property
    def h_vec(self):
        """Specific angular momentum vector. """
        h_vec = np.cross(self.r.to(u.km).value, self.v.to(u.km / u.s)) * u.km ** 2 / u.s
        return h_vec

    @property
    def arglat(self):
        """Argument of latitude. """
        arglat = (self.argp + self.nu) % (360 * u.deg)
        return arglat

    def rv(self):
        """Position and velocity vectors. """
        return self.r, self.v

    def coe(self):
        """Classical orbital elements. """
        return (
            self.a,
            self.ecc,
            self.inc.to(u.deg),
            self.raan.to(u.deg),
            self.argp.to(u.deg),
            self.nu.to(u.deg),
        )

    def pqw(self):
        """Perifocal frame (PQW) vectors. """
        if self.ecc < 1e-8:
            if abs(self.inc.to(u.rad).value) > 1e-8:
                node = np.cross([0, 0, 1], self.h_vec) / norm(self.h_vec)
                p_vec = node / norm(node)  # Circular inclined
            else:
                p_vec = [1, 0, 0] * u.one  # Circular equatorial
        else:
            p_vec = self.e_vec / self.ecc
        w_vec = self.h_vec / norm(self.h_vec)
        q_vec = np.cross(w_vec, p_vec) * u.one
        return p_vec, q_vec, w_vec

    def to_vectors(self):
        """Converts to position and velocity vector representation.

        """
        raise NotImplementedError

    def to_classical(self):
        """Converts to classical orbital elements representation.

        """
        raise NotImplementedError

    def to_equinoctial(self):
        """Converts to modified equinoctial elements representation.

        """
        raise NotImplementedError


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

        return RVState(self.attractor, r * u.km, v * u.km / u.s)

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

        return ModifiedEquinoctialState(
            self.attractor,
            p * u.km,
            f * u.rad,
            g * u.rad,
            h * u.rad,
            k * u.rad,
            L * u.rad,
        )


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
        (p, ecc, inc, raan, argp, nu) = rv2coe(
            self.attractor.k.to(u.km ** 3 / u.s ** 2).value,
            self.r.to(u.km).value,
            self.v.to(u.km / u.s).value,
        )

        return ClassicalState(
            self.attractor,
            p * u.km,
            ecc * u.one,
            inc * u.rad,
            raan * u.rad,
            argp * u.rad,
            nu * u.rad,
        )


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

        return ClassicalState(
            self.attractor,
            p * u.km,
            ecc * u.one,
            inc * u.rad,
            raan * u.rad,
            argp * u.rad,
            nu * u.rad,
        )
