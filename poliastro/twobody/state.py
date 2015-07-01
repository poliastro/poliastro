# coding: utf-8
"""Two body problem.

TODO
----
* Include other element sets and non-elliptic orbits.
* Improve consistency of units.

"""

import numpy as np

from astropy import time
from astropy import units as u

from poliastro.twobody.propagation import kepler
from poliastro.twobody.conversion import coe2rv, rv2coe

from poliastro.util import check_units, norm

J2000 = time.Time("J2000", scale='utc')


class State(object):
    """Class to represent the position of a body wrt to an attractor.

    """
    def __init__(self, attractor, epoch):
        """Constructor. To create a `State` object better use `from_vectors`
        and `from_classical` methods.

        Parameters
        ----------
        attractor : Body
            Main attractor.
        epoch : Time
            Epoch.

        """
        self.attractor = attractor
        self.epoch = epoch

    @staticmethod
    def from_vectors(attractor, r, v, epoch=J2000):
        """Return `State` object from position and velocity vectors.

        Parameters
        ----------
        attractor : Body
            Main attractor.
        r : Quantity
            Position vector wrt attractor center.
        v : Quantity
            Velocity vector.
        epoch : Time, optional
            Epoch, default to J2000.

        """
        if not check_units((r, v), (u.m, u.m / u.s)):
            raise u.UnitsError("Units must be consistent")

        return _RVState(attractor, r, v, epoch)

    @staticmethod
    def from_classical(attractor, a, ecc, inc, raan, argp, nu,
                       epoch=J2000):
        """Return `State` object from orbital elements.

        Parameters
        ----------
        attractor : Body
            Main attractor.
        a : Quantity
            Semimajor axis.
        ecc : Quantity
            Eccentricity.
        inc : Quantity
            Inclination
        raan : Quantity
            Right ascension of the ascending node.
        argp : Quantity
            Argument of the pericenter.
        nu : Quantity
            True anomaly.
        epoch : Time, optional
            Epoch, default to J2000.

        """
        if not check_units((a, ecc, inc, raan, argp, nu),
                           (u.m, u.one, u.rad, u.rad, u.rad, u.rad)):
            raise u.UnitsError("Units must be consistent")

        if ecc == 1.0 * u.one:
            raise ValueError("For parabolic orbits use "
                             "State.parabolic instead")

        ss = _ClassicalState(attractor, a, ecc, inc, raan, argp, nu, epoch)
        return ss

    @classmethod
    def circular(cls, attractor, alt,
                 inc=0 * u.deg, raan=0 * u.deg, arglat=0 * u.deg, epoch=J2000):
        """Return `State` corresponding to a circular orbit.

        Parameters
        ----------
        attractor : Body
            Main attractor.
        alt : Quantity
            Altitude over surface.
        inc : Quantity, optional
            Inclination, default to 0 deg (equatorial orbit).
        raan : Quantity, optional
            Right ascension of the ascending node, default to 0 deg.
        arglat : Quantity, optional
            Argument of latitude, default to 0 deg.
        epoch: Time, optional
            Epoch, default to J2000.

        """
        if not check_units((alt, inc, raan, arglat),
                           (u.m, u.rad, u.rad, u.rad)):
            raise u.UnitsError("Units must be consistent")
        a = attractor.R + alt
        ecc = 0 * u.one
        argp = 0 * u.deg
        ss = cls.from_classical(attractor, a, ecc, inc, raan, argp, arglat,
                                epoch)
        return ss

    @classmethod
    def parabolic(cls, attractor, p, inc, raan, argp, nu, epoch=J2000):
        """Return `State` corresponding to parabolic orbit.

        Parameters
        ----------
        attractor : Body
            Main attractor.
        p : Quantity
            Semi-latus rectum or parameter.
        inc : Quantity, optional
            Inclination.
        raan : Quantity
            Right ascension of the ascending node.
        argp : Quantity
            Argument of the pericenter.
        nu : Quantity
            True anomaly.
        epoch: Time, optional
            Epoch, default to J2000.

        """
        if not check_units((p, inc, raan, argp, nu),
                           (u.m, u.rad, u.rad, u.rad, u.rad)):
            raise u.UnitsError("Units must be consistent")

        k = attractor.k.to(u.km ** 3 / u.s ** 2)
        ecc = 1.0 * u.one
        r, v = coe2rv(k.to(u.km ** 3 / u.s ** 2).value,
                      p.to(u.km).value, ecc.value, inc.to(u.rad).value,
                      raan.to(u.rad).value, argp.to(u.rad).value,
                      nu.to(u.rad).value)

        ss = cls.from_vectors(attractor, r * u.km, v * u.km / u.s, epoch)
        return ss

    @property
    def r(self):
        """Position vector. """
        return self.to_vectors().r

    @property
    def v(self):
        """Velocity vector. """
        return self.to_vectors().v

    @property
    def coe(self):
        """Classical orbital elements. """
        return self.a, self.ecc, self.inc, self.raan, self.argp, self.nu

    @property
    def a(self):
        """Semimajor axis. """
        return self.to_classical().a

    @property
    def p(self):
        """Semilatus rectum. """
        return self.a * (1 - self.ecc ** 2)

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
    def period(self):
        """Period of the orbit. """
        return 2 * np.pi * u.rad / self.n

    @property
    def n(self):
        """Mean motion. """
        return np.sqrt(self.attractor.k / self.a ** 3) * u.rad

    @property
    def h(self):
        """Specific angular momentum. """
        return norm(self.h_vec)

    @property
    def e_vec(self):
        """Eccentricity vector. """
        r, v = self.rv()
        k = self.attractor.k
        e_vec = ((v.dot(v) - k / (norm(r))) * r - r.dot(v) * v) / k
        return e_vec

    @property
    def h_vec(self):
        """Specific angular momentum vector. """
        h_vec = np.cross(self.r.to(u.km).value,
                         self.v.to(u.km / u.s)) * u.km ** 2 / u.s
        return h_vec

    def rv(self):
        """Position and velocity vectors. """
        return self.r, self.v

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
        w_vec = self.h_vec / self.h
        q_vec = np.cross(w_vec, p_vec) * u.one
        return p_vec, q_vec, w_vec

    def propagate(self, time_of_flight):
        """Propagate this `State` some `time` and return the result.

        """
        r, v = kepler(self.attractor.k.to(u.km ** 3 / u.s ** 2).value,
                      self.r.to(u.km).value, self.v.to(u.km / u.s).value,
                      time_of_flight.to(u.s).value)
        return self.from_vectors(self.attractor, r * u.km, v * u.km / u.s,
                                 self.epoch + time_of_flight)

    def apply_maneuver(self, maneuver, intermediate=False):
        """Returns resulting State after applying maneuver to self.

        Optionally return intermediate states (default to False).

        Parameters
        ----------
        maneuver : Maneuver
            Maneuver to apply.
        intermediate : bool, optional
            Return intermediate states, default to False.

        """
        ss_new = self  # Initialize
        states = []
        attractor = self.attractor
        for delta_t, delta_v in maneuver:
            if not delta_t == 0 * u.s:
                ss_new = ss_new.propagate(time_of_flight=delta_t)
            r, v = ss_new.rv()
            vnew = v + delta_v
            ss_new = ss_new.from_vectors(attractor, r, vnew, ss_new.epoch)
            states.append(ss_new)
        if intermediate:
            res = states
        else:
            res = ss_new
        return res


class _RVState(State):
    def __init__(self, attractor, r, v, epoch):
        super(_RVState, self).__init__(attractor, epoch)
        self._r = r
        self._v = v

    @property
    def r(self):
        return self._r

    @property
    def v(self):
        return self._v

    def to_vectors(self):
        return self

    def to_classical(self):
        (a, ecc, inc, raan, argp, nu
         ) = rv2coe(self.attractor.k.to(u.km ** 3 / u.s ** 2).value,
                    self.r.to(u.km).value,
                    self.v.to(u.km / u.s).value)
        return super(_RVState, self).from_classical(self.attractor,
                                                    a * u.km,
                                                    ecc * u.one,
                                                    inc * u.rad,
                                                    raan * u.rad,
                                                    argp * u.rad,
                                                    nu * u.rad,
                                                    self.epoch)


class _ClassicalState(State):
    def __init__(self, attractor, a, ecc, inc, raan, argp, nu,
                 epoch):
        super(_ClassicalState, self).__init__(attractor, epoch)
        self._a = a
        self._ecc = ecc
        self._inc = inc
        self._raan = raan
        self._argp = argp
        self._nu = nu

    @property
    def a(self):
        return self._a

    @property
    def ecc(self):
        return self._ecc

    @property
    def inc(self):
        return self._inc

    @property
    def raan(self):
        return self._raan

    @property
    def argp(self):
        return self._argp

    @property
    def nu(self):
        return self._nu

    def _repr_latex_(self):
        """Creates a LaTeX representation.

        Used by the IPython notebook.

        """
        elem_names = [r"a", r"e", r"i", r"\Omega", r"\omega", r"\nu"]
        elem_values = [elem._repr_latex_().strip("$")
                       for elem in self.elements]
        pairs = zip(elem_names, elem_values)
        res = r"\\".join(["{0} & = {1}".format(name, value)
                         for name, value in pairs])
        return r"$\begin{{align}}{}\end{{align}}$".format(res)

    def to_vectors(self):
        r, v = coe2rv(self.attractor.k.to(u.km ** 3 / u.s ** 2).value,
                      self.p.to(u.km).value,
                      self.ecc.value,
                      self.inc.to(u.rad).value,
                      self.raan.to(u.rad).value,
                      self.argp.to(u.rad).value,
                      self.nu.to(u.rad).value)

        return super(_ClassicalState, self).from_vectors(self.attractor,
                                                         r * u.km,
                                                         v * u.km / u.s,
                                                         self.epoch)

    def to_classical(self):
        return self
