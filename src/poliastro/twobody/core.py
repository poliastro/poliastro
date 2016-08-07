# coding: utf-8
"""Two body problem.

TODO
----
* Improve consistency of units.

"""
import numpy as np

from astropy import time
from astropy import units as u

from poliastro.twobody.propagation import kepler

from poliastro.util import norm

J2000 = time.Time("J2000", scale='utc')


class State(object):
    """Class to represent the position of a body wrt to an attractor.

    """
    def __init__(self, attractor, epoch):
        """Constructor. To create a `State` object use `from_vectors`
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
    @u.quantity_input(r=u.m, v=u.m / u.s)
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
        assert np.any(r.value), "Position vector must be non zero"

        return poliastro.twobody.rv.RVState(
            attractor, r, v, epoch)

    @staticmethod
    @u.quantity_input(a=u.m, ecc=u.one, inc=u.rad, raan=u.rad, argp=u.rad, nu=u.rad)
    def from_classical(attractor, a, ecc, inc, raan, argp, nu,
                       epoch=J2000):
        """Return `State` object from classical orbital elements.

        Parameters
        ----------
        attractor : Body
            Main attractor.
        p : Quantity
            Semilatus rectum.
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
        if ecc == 1.0 * u.one:
            raise ValueError("For parabolic orbits use "
                             "State.parabolic instead")

        return poliastro.twobody.classical.ClassicalState(
            attractor, a, ecc, inc, raan, argp, nu, epoch)

    @staticmethod
    @u.quantity_input(p=u.m, f=u.one, g=u.rad, h=u.rad, k=u.rad, L=u.rad)
    def from_equinoctial(attractor, p, f, g, h, k, L, epoch=J2000):
        """Return `State` object from modified equinoctial elements.

        Parameters
        ----------
        attractor : Body
            Main attractor.
        p : Quantity
            Semilatus rectum.
        f : Quantity
            Second modified equinoctial element.
        g : Quantity
            Third modified equinoctial element.
        h : Quantity
            Fourth modified equinoctial element.
        k : Quantity
            Fifth modified equinoctial element.
        L : Quantity
            True longitude.
        epoch : Time, optional
            Epoch, default to J2000.

        """
        return poliastro.twobody.equinoctial.ModifiedEquinoctialState(
            attractor, p, f, g, h, k, L, epoch)

    @classmethod
    @u.quantity_input(alt=u.m, inc=u.rad, raan=u.rad, arglat=u.rad)
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
        a = attractor.R + alt
        ecc = 0 * u.one
        argp = 0 * u.deg

        return cls.from_classical(attractor, a, ecc, inc, raan, argp, arglat,
                                  epoch)

    @classmethod
    @u.quantity_input(p=u.m, inc=u.rad, raan=u.rad, argp=u.rad, nu=u.rad)
    def parabolic(cls, attractor, p, inc, raan, argp, nu, epoch=J2000):
        """Return `State` corresponding to parabolic orbit.

        Parameters
        ----------
        attractor : Body
            Main attractor.
        p : Quantity
            Semilatus rectum or parameter.
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
        k = attractor.k.to(u.km ** 3 / u.s ** 2)
        ecc = 1.0 * u.one
        r, v = poliastro.twobody.classical.coe2rv(
                k.to(u.km ** 3 / u.s ** 2).value,
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
    def a(self):
        """Semimajor axis. """
        return self.to_classical().a

    @property
    def p(self):
        """Semilatus rectum. """
        return self.a * (1 - self.ecc**2)

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
        return np.sqrt(self.attractor.k / self.a ** 3) * u.rad

    @property
    def energy(self):
        """Specific energy. """
        return (self.v.dot(self.v) / 2 -
                self.attractor.k / np.sqrt(self.r.dot(self.r)))

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
        h_vec = np.cross(self.r.to(u.km).value,
                         self.v.to(u.km / u.s)) * u.km ** 2 / u.s
        return h_vec

    @property
    def arglat(self):
        """Argument of latitude. """
        arglat = self.argp + self.nu
        return arglat

    def rv(self):
        """Position and velocity vectors. """
        return self.r, self.v

    def coe(self):
        """Classical orbital elements. """
        return self.a, self.ecc, self.inc, self.raan, self.argp, self.nu

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

    def propagate(self, time_of_flight, rtol=1e-10):
        """Propagate this `State` some `time` and return the result.

        """
        r, v = kepler(self.attractor.k.to(u.km ** 3 / u.s ** 2).value,
                      self.r.to(u.km).value, self.v.to(u.km / u.s).value,
                      time_of_flight.to(u.s).value,
                      rtol=rtol)
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


# Imports at the bottom to avoid circular import problems
import poliastro.twobody.rv
import poliastro.twobody.classical
import poliastro.twobody.equinoctial
