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
from astropy import units as u
from poliastro import constants


class _Body(object):
    def __str__(self):
        return u"{0} ({1})".format(self.name, self.symbol)

    def __repr__(self):
        return self.__str__()


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


class _Venus(_Body):
    parent = Sun
    k = constants.GM_venus
    name = "Venus"
    symbol = u"\u2640"
    R = constants.R_venus


class _Earth(_Body):
    parent = Sun
    k = constants.GM_earth
    name = "Earth"
    symbol = u"\u2641"
    R = constants.R_earth


class _Mars(_Body):
    parent = Sun
    k = constants.GM_mars
    name = "Mars"
    symbol = u"\u2642"
    R = constants.R_mars


class _Jupiter(_Body):
    parent = Sun
    k = constants.GM_jupiter
    name = "Jupiter"
    symbol = u"\u2643"
    R = constants.R_jupiter


class _Saturn(_Body):
    parent = Sun
    k = constants.GM_saturn
    name = "Saturn"
    symbol = u"\u2644"
    R = constants.R_saturn


class _Uranus(_Body):
    parent = Sun
    k = constants.GM_uranus
    name = "Uranus"
    symbol = u"\u26E2"
    R = constants.R_uranus


class _Neptune(_Body):
    parent = Sun
    k = constants.GM_neptune
    name = "Neptune"
    symbol = u"\u2646"
    R = constants.R_neptune


class _Pluto(_Body):
    parent = Sun
    k = constants.GM_pluto
    name = "Pluto"
    symbol = u"\u2647"
    R = constants.R_pluto

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
