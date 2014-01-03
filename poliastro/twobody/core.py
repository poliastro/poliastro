# coding: utf-8
"""Two body problem.

"""

import numpy as np
from numpy.linalg import norm

from astropy import time
from astropy import units as u

from poliastro.util import transform

J2000 = time.Time("J2000", scale='utc')


class State(object):
    """Class to represent a position of a body wrt to an attractor.

    """
    def __init__(self, attractor, values, epoch=J2000):
        """Constructor.

        Parameters
        ----------
        attractor : Body
            Main attractor.
        values : tuple
            Values to define the state, either orbital elements
            or r, v vectors.
        epoch : Time, optional
            Epoch, default to J2000

        """
        # TODO: Cannot define a parabola with its semi-major axis
        self.attractor = attractor
        self.epoch = epoch
        if len(values) == 6:
            _check_elements_units(values)
        elif len(values) == 2:
            _check_rv_units(values)
        else:
            raise ValueError("Incorrect number of parameters")
        self._values = values

    @classmethod
    def coe2rv(cls, attractor, elements):
        """Return `State` object using r and v from elements.

        """
        # TODO: Desirable?
        #ss_coe = cls(attractor, elements)
        #ss_coe.p, ss_coe.ecc, ...
        a, ecc, inc, raan, argp, nu = elements
        p = a * (1 - ecc ** 2)
        r_pqw = [np.cos(nu) / (1 + ecc * np.cos(nu)),
                 np.sin(nu) / (1 + ecc * np.cos(nu)),
                 0] * p
        v_pqw = [-np.sin(nu),
                 (ecc + np.cos(nu)),
                 0] * np.sqrt(attractor.k / p).to(u.km / u.s)
        r_ijk = transform(r_pqw, -argp, 'z')
        r_ijk = transform(r_ijk, -inc, 'x')
        r_ijk = transform(r_ijk, -raan, 'z')
        v_ijk = transform(v_pqw, -argp, 'z')
        v_ijk = transform(v_ijk, -inc, 'x')
        v_ijk = transform(v_ijk, -raan, 'z')

        return cls(attractor, (r_ijk, v_ijk))

    @classmethod
    def rv2coe(cls, attractor, rv):
        """Returns `State` object using elements from r and v.

        """
        # TODO: See coe2rv
        # HACK: Neither np.dot nor np.cross won't preserve Quantity,
        # see astropy #1930

        # Initial data
        r, v = rv
        k = attractor.k

        # Shape parameters
        h = np.cross(r.to(u.km), v.to(u.km / u.s)) * u.km ** 2 / u.s
        n = np.cross([0, 0, 1], h) / norm(h)
        e = ((((v.dot(v) - k / (norm(r) * r.unit)) * r - r.dot(v) * v) / k)
             .decompose().value)
        ecc = norm(e)
        p = h.dot(h) / k
        a = p / (1 - ecc ** 2)  # Might produce np.inf

        # Angular parameters
        inc = np.arccos(h[2] / (norm(h) * h.unit)).to(u.deg)
        raan = (np.arctan2(n[1], n[0])
                % (2 * np.pi) * u.rad).to(u.deg)
        argp = (np.arctan2(h.value.dot(np.cross(n, e)) / norm(h), e.dot(n))
                % (2 * np.pi) * u.rad).to(u.deg)
        nu = (np.arctan2(h.value.dot(np.cross(e, r)) / norm(h), r.value.dot(e))
              % (2 * np.pi) * u.rad).to(u.deg)

        return cls(attractor, (a, ecc, inc, raan, argp, nu))

    @property
    def elements(self):
        """Classical orbital elements.

        """
        if len(self._values) == 6:
            return self._values
        else:
            return self.rv2coe(self.attractor, self.rv).elements

    @property
    def rv(self):
        """Position and velocity vectors.

        """
        if len(self._values) == 2:
            return self._values
        else:
            return self.coe2rv(self.attractor, self.elements).rv


def _check_elements_units(elements):
    """Check if orbital elements have consistent units.

    """
    ELEMENTS_UNITS = (u.m, u.dimensionless_unscaled, u.rad,
                      u.rad, u.rad, u.rad)
    for ii, unit in enumerate(ELEMENTS_UNITS):
        try:
            if elements[ii].si.unit != unit:
                raise u.UnitsError("Units must be consistent")
        except AttributeError:
            if ii != 1:
                raise ValueError("Elements must have units "
                                 "(use astropy.units)")


def _check_rv_units(rv):
    """Check if r and v vectors have consistent units.

    """
    r, v = rv
    try:
        if r.si.unit != u.m or v.si.unit != u.m / u.s:
            raise u.UnitsError("Units must be consistent")
    except AttributeError:
        raise ValueError("r and v vectors must have units "
                         "(use astropy.units)")
