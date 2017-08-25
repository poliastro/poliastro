# coding: utf-8
"""Bodies of the Solar System.

Contains some predefined bodies of the Solar System:

* Sun (☉)
* Earth (♁)
* Moon (☾)
* Mercury (☿)
* Venus (♀)
* Mars (♂)
* Jupiter (♃)
* Saturn (♄)
* Uranus (⛢)
* Neptune (♆)
* Pluto (♇)

and a way to define new bodies (:py:class:`~Body` class).

Data references can be found in :py:mod:`~poliastro.constants`
"""
import math
from functools import wraps

from astropy import units as u
from poliastro import constants


def d_and_t_from_epoch(func):

    @wraps(func)
    def wrapper(epoch):
        T = (epoch.tdb - constants.J2000).to('day').value / 36525
        d = (epoch.tdb - constants.J2000).to('day').value

        return func(T, d)

    return wrapper


class _Body(object):
    def __str__(self):
        return u"{0} ({1})".format(self.name, self.symbol)

    def __repr__(self):
        return self.__str__()

    @staticmethod
    def rot_elements_at_epoch(epoch):
        """Provides rotational elements at epoch.

        Provides north pole of body and angle to prime meridian

        Parameters
        ----------
        epoch : ~astropy.time.Time, optional
            Epoch, default to J2000.

        Returns
        -------
        ra, dec, W: tuple (~astropy.units.Quantity)
            Right ascension and declination of north pole, and angle of the prime meridian.

        """
        raise NotImplementedError('Function only defined for some Solar System bodies')


class Body(_Body):
    """Class to represent a generic body of the Solar System.

    """
    def __init__(self, parent, k, name, symbol=None, R=0 * u.km):
        """Constructor.

        Parameters
        ----------
        parent : Body
            Central body.
        k : ~astropy.units.Quantity
            Standard gravitational parameter.
        name : str
            Name of the body.
        symbol : str, optional
            Symbol for the body.
        R : ~astropy.units.Quantity, optional
            Radius of the body.

        """
        self.parent = parent
        self.k = k
        self.name = name
        self.symbol = symbol
        self.R = R

    @classmethod
    @u.quantity_input(k=u.km ** 3 / u.s ** 2, R=u.km)
    def from_parameters(cls, parent, k, name, symbol, R):
        return cls(parent, k, name, symbol, R)


class _Sun(_Body):
    parent = None
    k = constants.GM_sun
    name = "Sun"
    symbol = u"\u2609"
    R = constants.R_sun

Sun = _Sun()


class _Mercury(_Body):
    parent = Sun
    k = constants.GM_mercury
    name = "Mercury"
    symbol = u"\u263F"
    R = constants.R_mercury

    @staticmethod
    @d_and_t_from_epoch
    def rot_elements_at_epoch(T, d):
        ra = (281.01 - 0.033 * T) * u.deg
        dec = (61.45 - 0.005 * T) * u.deg
        W = (329.548 + 6.1385025 * d) * u.deg

        return ra, dec, W


class _Venus(_Body):
    parent = Sun
    k = constants.GM_venus
    name = "Venus"
    symbol = u"\u2640"
    R = constants.R_venus

    @staticmethod
    @d_and_t_from_epoch
    def rot_elements_at_epoch(T, d):

        ra = 272.76 * u.deg
        dec = 67.16 * u.deg
        W = (160.20 - 1.4813688 * d) * u.deg

        return ra, dec, W


class _Earth(_Body):
    parent = Sun
    k = constants.GM_earth
    name = "Earth"
    symbol = u"\u2641"
    R = constants.R_earth

    @staticmethod
    @d_and_t_from_epoch
    def rot_elements_at_epoch(T, d):

        ra = (0.00 - 0.641 * T) * u.deg
        dec = (90.00 - 0.557 * T) * u.deg
        W = (190.147 + 360.9856235 * d) * u.deg

        return ra, dec, W


class _Mars(_Body):
    parent = Sun
    k = constants.GM_mars
    name = "Mars"
    symbol = u"\u2642"
    R = constants.R_mars

    @staticmethod
    @d_and_t_from_epoch
    def rot_elements_at_epoch(T, d):

        ra = (317.68143 - 0.1061 * T) * u.deg
        dec = (52.88650 - 0.0609 * T) * u.deg
        W = (176.630 + 350.89198226 * d) * u.deg

        return ra, dec, W


class _Jupiter(_Body):
    parent = Sun
    k = constants.GM_jupiter
    name = "Jupiter"
    symbol = u"\u2643"
    R = constants.R_jupiter

    @staticmethod
    @d_and_t_from_epoch
    def rot_elements_at_epoch(T, d):

        Ja = (99.360714 + 4850.4046 * T) * u.deg
        Jb = (175.895369 + 1191.9605 * T) * u.deg
        Jc = (300.323162 + 262.5475 * T) * u.deg
        Jd = (114.012305 + 6070.2476 * T) * u.deg
        Je = (49.511251 + 64.3000 * T) * u.deg

        ra = (268.056595 - 0.006499 * T + 0.000117 * math.sin(Ja) + 0.000938 * math.sin(Jb) +
              0.001432 * math.sin(Jc) + 0.000030 * math.sin(Jd) + 0.002150 * math.sin(Je)) * u.deg
        dec = (64.495303 + 0.002413 * T + 0.000050 * math.cos(Ja) + 0.000404 * math.cos(Jb) +
               0.000617 * math.cos(Jc) - 0.000013 * math.cos(Jd) + 0.000926 * math.cos(Je)) * u.deg
        W = (284.95 + 870.5366420 * d) * u.deg

        return ra, dec, W


class _Saturn(_Body):
    parent = Sun
    k = constants.GM_saturn
    name = "Saturn"
    symbol = u"\u2644"
    R = constants.R_saturn

    @staticmethod
    @d_and_t_from_epoch
    def rot_elements_at_epoch(T, d):

        ra = (40.589 - 0.036 * T) * u.deg
        dec = (83.537 - 0.004 * T) * u.deg
        W = (38.90 + 810.7939024 * d) * u.deg

        return ra, dec, W


class _Uranus(_Body):
    parent = Sun
    k = constants.GM_uranus
    name = "Uranus"
    symbol = u"\u26E2"
    R = constants.R_uranus

    @staticmethod
    @d_and_t_from_epoch
    def rot_elements_at_epoch(T, d):
        ra = 257.311 * u.deg
        dec = -15.175 * u.deg
        W = (203.81 - 501.1600928 * d) * u.deg

        return ra, dec, W


class _Neptune(_Body):
    parent = Sun
    k = constants.GM_neptune
    name = "Neptune"
    symbol = u"\u2646"
    R = constants.R_neptune

    @staticmethod
    @d_and_t_from_epoch
    def rot_elements_at_epoch(T, d):
        ra = (299.36 + 0.70 * math.sin(N)) * u.deg
        dec = (43.46 - 0.51 * math.cos(N)) * u.deg
        W = (253.18 + 536.3128492 * d - 0.48 * math.sin(N)) * u.deg

        return ra, dec, W


class _Pluto(_Body):
    parent = Sun
    k = constants.GM_pluto
    name = "Pluto"
    symbol = u"\u2647"
    R = constants.R_pluto

    @staticmethod
    @d_and_t_from_epoch
    def rot_elements_at_epoch(T, d):
        ra = 312.993 * u.deg
        dec = 6.163 * u.deg
        W = (190.147 + 360.9856235 * d) * u.deg

        return ra, dec, W


Mercury = _Mercury()
Venus = _Venus()
Earth = _Earth()
Mars = _Mars()
Jupiter = _Jupiter()
Saturn = _Saturn()
Uranus = _Uranus()
Neptune = _Neptune()

Pluto = _Pluto()


class _Moon(_Body):
    parent = Earth
    k = constants.GM_moon
    name = "Moon"
    symbol = u"\u263E"
    R = constants.R_moon

Moon = _Moon()
