# coding: utf-8
"""Functions to define orbits from classical orbital elements.

"""
import numpy as np
from numpy import sin, cos, sqrt
from astropy import units as u

from poliastro.util import transform

from poliastro.twobody.core import State


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

    Notes
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


class ClassicalState(State):
    def __init__(self, attractor, p, ecc, inc, raan, argp, nu,
                 epoch):
        super(ClassicalState, self).__init__(attractor, epoch)
        self._p = p
        self._ecc = ecc
        self._inc = inc
        self._raan = raan
        self._argp = argp
        self._nu = nu

    @property
    def p(self):
        return self._p

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
        # noinspection PyProtectedMember
        # Each element is an astropy Quantity:
        # https://github.com/astropy/astropy/blob/v1.0.7/astropy/units/quantity.py#L935
        # '_repr_index' is an internal method for IPython use:
        # http://ipython.readthedocs.org/en/3.x/config/integrating.html#rich-display
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

        return super(ClassicalState, self).from_vectors(self.attractor,
                                                        r * u.km,
                                                        v * u.km / u.s,
                                                        self.epoch)

    def to_classical(self):
        return self

    def to_equinoctial(self):
        p, f, g, h, k, L = coe2mee(self.p.to(u.km).value,
                                   self.ecc.value,
                                   self.inc.to(u.rad).value,
                                   self.raan.to(u.rad).value,
                                   self.argp.to(u.rad).value,
                                   self.nu.to(u.rad).value)

        return super(ClassicalState, self).from_equinoctial(
            self.attractor,
            p * u.km,
            f * u.rad,
            g * u.rad,
            h * u.rad,
            k * u.rad,
            L * u.rad,
            self.epoch)
