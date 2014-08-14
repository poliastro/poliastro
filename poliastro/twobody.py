# coding: utf-8
"""Two body problem.

TODO
----
* Include other element sets and non-elliptic orbits.
* Improve consistency of units.

"""

import numpy as np
from numpy.linalg import norm

from astropy import time
from astropy import units as u
u.one = u.dimensionless_unscaled  # astropy #1980

from poliastro.util import transform, check_units

# For plotting
import matplotlib.pyplot as plt
import matplotlib as mpl

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

    def plot2D(self, ax=None, num=100, osculating=True):
        """Plots state and osculating orbit in their plane.

        Parameters
        ----------
        ax : Axes, optional
            Matplotlib axes where plotting should be performed. If not given,
            new ones will be created.
        num : int, optional
            Number of points, default to 100.
        osculating : bool, optional
            Whether to plot the osculating orbit, default to yes.

        Notes
        -----
        TODO: I still miss a way to plot several orbits in one plot. Probably
        it's time to create some plotting.py module, with a OrbitPlotter where
        I can hold things like the plane.

        """
        if not ax:
            _, ax = plt.subplots(figsize=(6, 6))

        lines = []
        nu_vals = np.linspace(0, 2 * np.pi, num) + self.nu.to(u.rad).value
        r_pqw, _ = rv_pqw(self.attractor.k.to(u.km ** 3 / u.s ** 2).value,
                          self.p.to(u.km).value, self.ecc.value, nu_vals)

        # Current position
        l, = ax.plot(r_pqw[0, 0], r_pqw[0, 1], 'o')
        lines.append(l)

        # Attractor
        attractor = mpl.patches.Circle((0, 0),
                                       self.attractor.R.to(u.km).value,
                                       lw=0, color='#204a87')  # Earth
        ax.add_patch(attractor)

        if osculating:
            l, = ax.plot(r_pqw[:, 0], r_pqw[:, 1], '--', color=l.get_color())
            lines.append(l)

        # Appearance
        ax.set_title(self.epoch.iso)
        ax.set_xlabel("$x$ (km)")
        ax.set_ylabel("$y$ (km)")
        ax.set_aspect(1)

        return lines

    def plot(self):
        """Shortcut to `plot2D`.

        """
        return self.plot2D()


def rv_pqw(k, p, ecc, nu):
    """Returns r and v vectors in perifocal frame.

    """
    r_pqw = (np.array([np.cos(nu), np.sin(nu), 0 * nu]) * p / (1 + ecc * np.cos(nu))).T
    v_pqw = (np.array([-np.sin(nu), (ecc + np.cos(nu)), 0]) * np.sqrt(k / p)).T
    return r_pqw, v_pqw


def coe2rv(k, a, ecc, inc, raan, argp, nu):
    """Converts from orbital elements to vectors.

    Parameters
    ----------
    k : float
        Standard gravitational parameter (km^3 / s^2).
    a : float
        Semi-major axis (km).
    ecc : float
        Eccentricity.
    inc : float
        Inclination (rad).
    omega : float
        Longitude of ascending node (rad).
    argp : float
        Argument of perigee (rad).
    nu : float
        True anomaly (rad).

    """
    p = a * (1 - ecc ** 2)
    r_pqw, v_pqw = rv_pqw(k, p, ecc, nu)

    r_ijk = transform(r_pqw, -argp, 'z', u.rad)
    r_ijk = transform(r_ijk, -inc, 'x', u.rad)
    r_ijk = transform(r_ijk, -raan, 'z', u.rad)
    v_ijk = transform(v_pqw, -argp, 'z', u.rad)
    v_ijk = transform(v_ijk, -inc, 'x', u.rad)
    v_ijk = transform(v_ijk, -raan, 'z', u.rad)

    return r_ijk, v_ijk


def rv2coe(k, r, v):
    """Converts from vectors to orbital elements.

    Parameters
    ----------
    k : float
        Standard gravitational parameter (km^3 / s^2).
    r : array
        Position vector (km).
    v : array
        Velocity vector (km / s).

    """
    h = np.cross(r, v)
    n = np.cross([0, 0, 1], h) / norm(h)
    e = ((v.dot(v) - k / (norm(r))) * r - r.dot(v) * v) / k
    ecc = norm(e)
    p = h.dot(h) / k
    a = p / (1 - ecc ** 2)

    inc = np.arccos(h[2] / norm(h))
    raan = np.arctan2(n[1], n[0]) % (2 * np.pi)
    argp = np.arctan2(h.dot(np.cross(n, e)) / norm(h), e.dot(n)) % (2 * np.pi)
    nu = np.arctan2(h.dot(np.cross(e, r)) / norm(h), r.dot(e)) % (2 * np.pi)

    return a, ecc, inc, raan, argp, nu


def kepler(k, r0, v0, tof):
    """Propagates orbit.

    This is a wrapper around kepler from ast2body.for.

    Parameters
    ----------
    k : float
        Gravitational constant of main attractor (km^3 / s^2).
    r0 : array
        Initial position (km).
    v0 : array
        Initial velocity (km).
    tof : float
        Time of flight (s).

    Raises
    ------
    RuntimeError
        If the status of the subroutine is not 'ok'.

    """
    r0 = np.asarray(r0).astype(np.float)
    v0 = np.asarray(v0).astype(np.float)
    tof = float(tof)
    assert r0.shape == (3,)
    assert v0.shape == (3,)
    r, v, error = _ast2body.kepler(r0, v0, tof, k)
    error = error.strip().decode('ascii')
    if error != 'ok':
        raise RuntimeError("There was an error: {}".format(error))
    return r, v
