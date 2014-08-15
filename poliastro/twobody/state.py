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
u.one = u.dimensionless_unscaled  # astropy #1980

from poliastro.twobody.propagation import kepler
from poliastro.twobody.conversion import coe2rv, rv2coe

from poliastro.plotting import OrbitPlotter

from poliastro.util import check_units

from . import _ast2body

J2000 = time.Time("J2000", scale='utc')


class State(object):
    """Class to represent the position of a body wrt to an attractor.

    """
    def __init__(self, attractor, r, v, epoch):
        """Constructor. To create a `State` object better use `from_vectors`
        and `from_elements` methods.

        Parameters
        ----------
        attractor : Body
            Main attractor.
        r, v : array
            Position and velocity vectors.
        epoch : Time
            Epoch.

        """
        self.attractor = attractor
        self.epoch = epoch
        self.r = r
        self.v = v
        self._elements = None

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

    @classmethod
    def from_vectors(cls, attractor, r, v, epoch=J2000):
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

        return cls(attractor, r, v, epoch)

    @classmethod
    def from_elements(cls, attractor, a, ecc, inc, raan, argp, nu,
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

        k = attractor.k.to(u.km ** 3 / u.s ** 2)
        r, v = coe2rv(k.to(u.km ** 3 / u.s ** 2).value,
                      a.to(u.km).value, ecc.value, inc.to(u.rad).value,
                      raan.to(u.rad).value, argp.to(u.rad).value,
                      nu.to(u.rad).value)

        ss = cls(attractor, r * u.km, v * u.km / u.s, epoch)
        ss._elements = a, ecc, inc, raan, argp, nu
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
        ss = cls.from_elements(attractor, a, ecc, inc, raan, argp, arglat,
                               epoch)
        return ss

    @property
    def elements(self):
        """Classical orbital elements.

        """
        if self._elements:
            return self._elements
        else:
            k = self.attractor.k.to(u.km ** 3 / u.s ** 2).value
            r = self.r.to(u.km).value
            v = self.v.to(u.km / u.s).value
            a, ecc, inc, raan, argp, nu = rv2coe(k, r, v)
            self._elements = (a * u.km, ecc * u.one, (inc * u.rad).to(u.deg),
                              (raan * u.rad).to(u.deg),
                              (argp * u.rad).to(u.deg),
                              (nu * u.rad).to(u.deg))
            return self._elements

    @property
    def a(self):
        """Semimajor axis.

        """
        a = self.elements[0]
        return a

    @property
    def p(self):
        """Semilatus rectum.

        """
        p = self.a * (1 - self.ecc ** 2)
        return p

    @property
    def r_p(self):
        """Radius of pericenter.

        """
        r_p = self.a * (1 - self.ecc)
        return r_p

    @property
    def r_a(self):
        """Radius of apocenter.

        """
        r_a = self.a * (1 + self.ecc)
        return r_a

    @property
    def ecc(self):
        """Eccentricity.

        """
        ecc = self.elements[1]
        return ecc

    @property
    def inc(self):
        """Inclination.

        """
        inc = self.elements[2]
        return inc

    @property
    def raan(self):
        """Right ascension of the ascending node.

        """
        raan = self.elements[3]
        return raan

    @property
    def argp(self):
        """Argument of the perigee.

        """
        argp = self.elements[4]
        return argp

    @property
    def nu(self):
        """True anomaly.

        """
        nu = self.elements[5]
        return nu

    @property
    def period(self):
        """Period of the orbit.

        """
        period = 2 * np.pi * u.rad / self.n
        return period

    @property
    def n(self):
        """Mean motion.

        """
        n = np.sqrt(self.attractor.k / self.a ** 3) * u.rad
        return n

    def rv(self):
        """Position and velocity vectors.

        """
        return self.r, self.v

    def propagate(self, time_of_flight):
        """Propagate this `State` some `time` and return the result.

        """
        r, v = kepler(self.attractor.k.to(u.km ** 3 / u.s ** 2).value,
                      self.r.to(u.km).value, self.v.to(u.km / u.s).value,
                      time_of_flight.to(u.s).value)
        return self.from_vectors(self.attractor, r * u.km, v * u.km / u.s,
                                 self.epoch + time_of_flight)

    def apply_maneuver(self, maneuver):
        """Returns resulting State after applying maneuver to self.

        Parameters
        ----------
        maneuver : Maneuver
            Maneuver to apply.

        """
        ss_new = self  # Initialize
        attractor = self.attractor
        for delta_t, delta_v in maneuver.impulses:
            if not delta_t == 0 * u.s:
                ss_new = ss_new.propagate(time_of_flight=delta_t)
            r, v = ss_new.rv()
            vnew = v + delta_v
            ss_new = ss_new.from_vectors(attractor, r, vnew)
        return ss_new

    def plot2D(self, ax=None, num=100, osculating=True):
        """Plots state and osculating orbit in their plane.

        This is a convenience function using
        :py:class:`poliastro.plotting.OrbitPlotter`.

        """
        op = OrbitPlotter(num)
        return op.plot(self, ax, osculating)

    def plot(self):
        """Shortcut to `plot2D`.

        """
        return self.plot2D()
