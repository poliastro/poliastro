# coding: utf-8
"""Two body problem.

TODO
----
* Improve consistency of units.

"""
import numpy as np

from astropy import time
from astropy import units as u

import poliastro.twobody.rv
import poliastro.twobody.classical
import poliastro.twobody.equinoctial

J2000 = time.Time("J2000", scale='utc')


class StateFactory(object):
    """Class to represent the position of a body wrt to an attractor.

    """
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
