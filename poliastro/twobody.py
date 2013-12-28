# coding: utf-8
"""Two body problem.

"""

import numpy as np

from astropy import units as u
from astropy import time

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
        # HACK: No easy way to build vector quantities, see
        # http://mail.scipy.org/pipermail/astropy/2013-December/002991.html
        r_pqw = [p.to(u.km).value * np.cos(nu) / (1 + ecc * np.cos(nu)),
                 p.to(u.km).value * np.sin(nu) / (1 + ecc * np.cos(nu)),
                 0] * u.km
        sqrtkp = np.sqrt(attractor.k / p).to(u.km / u.s)
        v_pqw = [-sqrtkp.value * np.sin(nu),
                 sqrtkp.value * (ecc + np.cos(nu)),
                 0] * u.km / u.s
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
        #import pdb; pdb.set_trace()
        # HACK: Neither np.dot nor np.cross won't preserve Quantity,
        # see astropy #1930
        def dot(a, b):
            return np.sum(a * b)
        r, v = rv
        k = attractor.k
        r_mag = np.sqrt(dot(r, r))
        v2 = dot(v, v)
        h = np.cross(r, v) * u.km ** 2 / u.s
        n = np.cross([0, 0, 1], h) * h.unit
        n_mag = np.sqrt(dot(n, n))
        e = (((v2 - k / r_mag) * r - dot(r, v) * v) / k).decompose()
        ecc = np.sqrt(dot(e, e))
        p = dot(h, h) / k
        a = p / (1 - ecc ** 2)  # Might produce np.inf
        inc = np.arccos(dot(h, [0, 0, 1]) / np.sqrt(dot(h, h))).to(u.deg)
        raan = np.arccos(dot(n, [1, 0, 0]) / np.sqrt(dot(n, n))).to(u.deg)
        # TODO: Use of arctan2 to automatically resolve quadrants
        if dot(n, [0, 1, 0]) < 0:
            raan = 360 * u.deg - raan
        argp = np.arccos(dot(n, e) / (ecc * n_mag)).to(u.deg)
        if dot(e, [0, 0, 1]) < 0:
            argp = 360 * u.deg - argp
        nu = np.arccos(dot(e, r) / (ecc * r_mag)).to(u.deg)
        if dot(r, v) < 0:
            nu = 360 * u.deg - nu
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
