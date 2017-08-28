# coding: utf-8
from datetime import datetime

import numpy as np

from astropy import units as u

from astropy import time

from poliastro.constants import J2000
from poliastro.ephem import get_body_ephem
from poliastro.twobody.propagation import propagate

import poliastro.twobody.rv
import poliastro.twobody.classical
import poliastro.twobody.equinoctial


ORBIT_FORMAT = "{r_p:.0f} x {r_a:.0f} x {inc:.1f} orbit around {body}"


class Orbit(object):
    """Position and velocity of a body with respect to an attractor
    at a given time (epoch).

    """
    def __init__(self, state, epoch):
        """Constructor.

        """
        #: Position and velocity or classical elements
        self.state = state
        #: Epoch of the orbit
        self.epoch = epoch

    @classmethod
    @u.quantity_input(r=u.m, v=u.m / u.s)
    def from_vectors(cls, attractor, r, v, epoch=J2000):
        """Return `Orbit` from position and velocity vectors.

        Parameters
        ----------
        attractor : Body
            Main attractor.
        r : ~astropy.units.Quantity
            Position vector wrt attractor center.
        v : ~astropy.units.Quantity
            Velocity vector.
        epoch : ~astropy.time.Time, optional
            Epoch, default to J2000.

        """
        assert np.any(r.value), "Position vector must be non zero"

        ss = poliastro.twobody.rv.RVState(
            attractor, r, v)
        return cls(ss, epoch)

    @classmethod
    @u.quantity_input(a=u.m, ecc=u.one, inc=u.rad, raan=u.rad, argp=u.rad, nu=u.rad)
    def from_classical(cls, attractor, a, ecc, inc, raan, argp, nu, epoch=J2000):
        """Return `Orbit` from classical orbital elements.

        Parameters
        ----------
        attractor : Body
            Main attractor.
        a : ~astropy.units.Quantity
            Semi-major axis.
        ecc : ~astropy.units.Quantity
            Eccentricity.
        inc : ~astropy.units.Quantity
            Inclination
        raan : ~astropy.units.Quantity
            Right ascension of the ascending node.
        argp : ~astropy.units.Quantity
            Argument of the pericenter.
        nu : ~astropy.units.Quantity
            True anomaly.
        epoch : ~astropy.time.Time, optional
            Epoch, default to J2000.

        """
        if ecc == 1.0 * u.one:
            raise ValueError("For parabolic orbits use Orbit.parabolic instead")
        elif not 0 * u.deg <= inc <= 180 * u.deg:
            raise ValueError("Inclination must be between 0 and 180 degrees")

        ss = poliastro.twobody.classical.ClassicalState(
            attractor, a, ecc, inc, raan, argp, nu)
        return cls(ss, epoch)

    @classmethod
    @u.quantity_input(p=u.m, f=u.one, g=u.rad, h=u.rad, k=u.rad, L=u.rad)
    def from_equinoctial(cls, attractor, p, f, g, h, k, L, epoch=J2000):
        """Return `Orbit` from modified equinoctial elements.

        Parameters
        ----------
        attractor : Body
            Main attractor.
        p : ~astropy.units.Quantity
            Semilatus rectum.
        f : ~astropy.units.Quantity
            Second modified equinoctial element.
        g : ~astropy.units.Quantity
            Third modified equinoctial element.
        h : ~astropy.units.Quantity
            Fourth modified equinoctial element.
        k : ~astropy.units.Quantity
            Fifth modified equinoctial element.
        L : ~astropy.units.Quantity
            True longitude.
        epoch : ~astropy.time.Time, optional
            Epoch, default to J2000.

        """
        ss = poliastro.twobody.equinoctial.ModifiedEquinoctialState(
            attractor, p, f, g, h, k, L)
        return cls(ss, epoch)

    @classmethod
    def from_body_ephem(cls, body, epoch=None):
        """Return osculating `Orbit` of a body at a given time.

        """
        if not epoch:
            epoch = time.Time.now()

        r, v = get_body_ephem(body.name, epoch)
        return cls.from_vectors(body.parent, r, v, epoch)

    @classmethod
    @u.quantity_input(alt=u.m, inc=u.rad, raan=u.rad, arglat=u.rad)
    def circular(cls, attractor, alt,
                 inc=0 * u.deg, raan=0 * u.deg, arglat=0 * u.deg, epoch=J2000):
        """Return circular `Orbit`.

        Parameters
        ----------
        attractor : Body
            Main attractor.
        alt : ~astropy.units.Quantity
            Altitude over surface.
        inc : ~astropy.units.Quantity, optional
            Inclination, default to 0 deg (equatorial orbit).
        raan : ~astropy.units.Quantity, optional
            Right ascension of the ascending node, default to 0 deg.
        arglat : ~astropy.units.Quantity, optional
            Argument of latitude, default to 0 deg.
        epoch: ~astropy.time.Time, optional
            Epoch, default to J2000.

        """
        a = attractor.R + alt
        ecc = 0 * u.one
        argp = 0 * u.deg

        return cls.from_classical(attractor, a, ecc, inc, raan, argp, arglat, epoch)

    @classmethod
    @u.quantity_input(p=u.m, inc=u.rad, raan=u.rad, argp=u.rad, nu=u.rad)
    def parabolic(cls, attractor, p, inc, raan, argp, nu, epoch=J2000):
        """Return parabolic `Orbit`.

        Parameters
        ----------
        attractor : Body
            Main attractor.
        p : ~astropy.units.Quantity
            Semilatus rectum or parameter.
        inc : ~astropy.units.Quantity, optional
            Inclination.
        raan : ~astropy.units.Quantity
            Right ascension of the ascending node.
        argp : ~astropy.units.Quantity
            Argument of the pericenter.
        nu : ~astropy.units.Quantity
            True anomaly.
        epoch: ~astropy.time.Time, optional
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

    def __str__(self):
        if self.a > 1e7 * u.km:
            unit = u.au
        else:
            unit = u.km

        return ORBIT_FORMAT.format(
            r_p=self.r_p.to(unit).value, r_a=self.r_a.to(unit), inc=self.inc.to(u.deg),
            body=self.attractor
        )

    def __repr__(self):
        return self.__str__()

    def propagate(self, time_of_flight, rtol=1e-10):
        """Propagate this `Orbit` some `time` and return the result.

        """
        return propagate(self, time_of_flight, rtol=rtol)

    def apply_maneuver(self, maneuver, intermediate=False):
        """Returns resulting `Orbit` after applying maneuver to self.

        Optionally return intermediate states (default to False).

        Parameters
        ----------
        maneuver : Maneuver
            Maneuver to apply.
        intermediate : bool, optional
            Return intermediate states, default to False.

        """
        orbit_new = self  # Initialize
        states = []
        attractor = self.attractor
        for delta_t, delta_v in maneuver:
            if not delta_t == 0 * u.s:
                orbit_new = orbit_new.propagate(time_of_flight=delta_t)
            r, v = orbit_new.rv()
            vnew = v + delta_v
            orbit_new = self.from_vectors(attractor, r, vnew, orbit_new.epoch)
            states.append(orbit_new)
        if intermediate:
            res = states
        else:
            res = orbit_new
        return res

    def rv(self):
        """Position and velocity vectors. """
        return self.state.rv()

    def coe(self):
        """Classical orbital elements. """
        return self.state.coe()

    def pqw(self):
        """Perifocal frame (PQW) vectors. """
        return self.state.pqw()

    @property
    def attractor(self):
        """Main attractor body. """
        return self.state.attractor

    @property
    def r(self):
        """Position vector. """
        return self.state.r

    @property
    def v(self):
        """Velocity vector. """
        return self.state.v

    @property
    def a(self):
        """Semimajor axis. """
        return self.state.a

    @property
    def p(self):
        """Semilatus rectum. """
        return self.state.p

    @property
    def r_p(self):
        """Radius of pericenter. """
        return self.state.r_p

    @property
    def r_a(self):
        """Radius of apocenter. """
        return self.state.r_a

    @property
    def ecc(self):
        """Eccentricity. """
        return self.state.ecc

    @property
    def inc(self):
        """Inclination. """
        return self.state.inc

    @property
    def raan(self):
        """Right ascension of the ascending node. """
        return self.state.raan

    @property
    def argp(self):
        """Argument of the perigee. """
        return self.state.argp

    @property
    def nu(self):
        """True anomaly. """
        return self.state.nu

    @property
    def f(self):
        """Second modified equinoctial element. """
        return self.state.f

    @property
    def g(self):
        """Third modified equinoctial element. """
        return self.state.g

    @property
    def h(self):
        """Fourth modified equinoctial element. """
        return self.state.h

    @property
    def k(self):
        """Fifth modified equinoctial element. """
        return self.state.k

    @property
    def L(self):
        """True longitude. """
        return self.state.L

    @property
    def period(self):
        """Period of the orbit. """
        return self.state.period

    @property
    def n(self):
        """Mean motion. """
        return self.state.n

    @property
    def energy(self):
        """Specific energy. """
        return self.state.energy

    @property
    def e_vec(self):
        """Eccentricity vector. """
        return self.state.e_vec

    @property
    def h_vec(self):
        """Specific angular momentum vector. """
        return self.state.h_vec

    @property
    def arglat(self):
        """Argument of latitude. """
        return self.state.arglat
