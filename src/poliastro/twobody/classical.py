"""Functions to define orbits from classical orbital elements.

"""
import numpy as np
from numpy import sin, cos, sqrt
from astropy import units as u

from poliastro.util import transform

from poliastro.twobody import rv
from poliastro.twobody import equinoctial

from ._base import BaseState


def rv_pqw(k, p, ecc, nu):
    """Returns r and v vectors in perifocal frame.

    """
    r_pqw = (np.array([cos(nu), sin(nu), 0 * nu]) * p / (1 + ecc * cos(nu))).T
    v_pqw = (np.array([-sin(nu), (ecc + cos(nu)), 0]) * sqrt(k / p)).T
    return r_pqw, v_pqw


def coe2rv(k, p, ecc, inc, raan, argp, nu):
    """Converts from classical orbital elements to vectors.

    Parameters
    ----------
    k : float
        Standard gravitational parameter (km^3 / s^2).
    p : float
        Semi-latus rectum or parameter (km).
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
    r_pqw, v_pqw = rv_pqw(k, p, ecc, nu)

    r_ijk = transform(r_pqw, -argp, 'z', u.rad)
    r_ijk = transform(r_ijk, -inc, 'x', u.rad)
    r_ijk = transform(r_ijk, -raan, 'z', u.rad)
    v_ijk = transform(v_pqw, -argp, 'z', u.rad)
    v_ijk = transform(v_ijk, -inc, 'x', u.rad)
    v_ijk = transform(v_ijk, -raan, 'z', u.rad)

    return r_ijk, v_ijk


def coe2mee(p, ecc, inc, raan, argp, nu):
    """Converts from classical orbital elements to modified equinoctial
    orbital elements.

    The definition of the modified equinoctial orbital elements is taken from
    [Walker, 1985].

    Parameters
    ----------
    k : float
        Standard gravitational parameter (km^3 / s^2).
    p : float
        Semi-latus rectum or parameter (km).
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

    Note
    -----
    The conversion equations are taken directly from the original paper.

    """
    lonper = raan + argp
    f = ecc * np.cos(lonper)
    g = ecc * np.sin(lonper)
    # TODO: Check polar case (see [Walker, 1985])
    h = np.tan(inc / 2) * np.cos(raan)
    k = np.tan(inc / 2) * np.sin(raan)
    L = lonper + nu
    return p, f, g, h, k, L


class ClassicalState(BaseState):
    """State defined by its classical orbital elements.

    """
    def __init__(self, attractor, p, ecc, inc, raan, argp, nu):
        super(ClassicalState, self).__init__(attractor)
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
        r, v = coe2rv(self.attractor.k.to(u.km ** 3 / u.s ** 2).value,
                      self.p.to(u.km).value,
                      self.ecc.value,
                      self.inc.to(u.rad).value,
                      self.raan.to(u.rad).value,
                      self.argp.to(u.rad).value,
                      self.nu.to(u.rad).value)

        return rv.RVState(
            self.attractor,
            r * u.km,
            v * u.km / u.s)

    def to_classical(self):
        """Converts to classical orbital elements representation.

        """
        return self

    def to_equinoctial(self):
        """Converts to modified equinoctial elements representation.

        """
        p, f, g, h, k, L = coe2mee(self.p.to(u.km).value,
                                   self.ecc.value,
                                   self.inc.to(u.rad).value,
                                   self.raan.to(u.rad).value,
                                   self.argp.to(u.rad).value,
                                   self.nu.to(u.rad).value)

        return equinoctial.ModifiedEquinoctialState(
            self.attractor,
            p * u.km,
            f * u.rad,
            g * u.rad,
            h * u.rad,
            k * u.rad,
            L * u.rad)
