# coding: utf-8
"""Bodies of the Solar System.

Contains some predefined bodies of the Solar System:

* Sun (☉)
* Earth (♁)

and a way to define new bodies (:py:class:`~Body` class).

"""
from astropy.constants import R_earth
from astropy import units as u


class Body(object):
    """Class to represent a body of the Solar System.

    """
    def __init__(self, parent, k, name, symbol=None, R=0 * u.km):
        """Constructor.

        Parameters
        ----------
        parent : Body
            Central body.
        k : Quantity
            Standard gravitational parameter.
        name : str
            Name of the body.
        symbol : str, optional
            Symbol for the body.
        R : Quantity, optional
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

    def __str__(self):
        return u"{0} ({1})".format(self.name, self.symbol)

    def _repr_latex_(self):
        """Creates a LaTeX representation.

        Used by the IPython notebook.

        """
        return self.__str__()


Sun = Body.from_parameters(
    None, k=132712440018 * u.km ** 3 / u.s ** 2,
    name="Sun", symbol=u"\u2609", R=695700 * u.km)
Mercury = Body.from_parameters(
    Sun, k=22032 * u.km ** 3 / u.s ** 2,
    name="Mercury", symbol=u"\u263F", R=2439.7 * u.km)
Venus = Body.from_parameters(
    Sun, k=324859 * u.km ** 3 / u.s ** 2,
    name="Venus", symbol=u"\u2640", R=6051.8 * u.km)
Earth = Body.from_parameters(
    Sun, k=398600 * u.km ** 3 / u.s ** 2,
    name="Earth", symbol=u"\u2641", R=R_earth.to(u.km))
Moon = Body.from_parameters(
    Earth, k=4904 * u.km ** 3 / u.s ** 2,
    name="Moon", symbol=u"\u1F315", R=1737.4 * u.km) # Used a white circle as symbol
Mars = Body.from_parameters(
    Sun, k=42828 * u.km ** 3 / u.s ** 2,
    name="Mars", symbol=u"\u2642", R=3389.5 * u.km)
Jupiter = Body.from_parameters(
    Sun, k=126686534 * u.km ** 3 / u.s ** 2,
    name="Jupiter", symbol=u"\u2643", R=69911 * u.km)
Saturn = Body.from_parameters(
    Sun, k=3793939 * u.km ** 3 / u.s ** 2,
    name="Saturn", symbol=u"\u2644", R=58232 * u.km)
Uranus = Body.from_parameters(
    Sun, k=5793939 * u.km ** 3 / u.s ** 2,
    name="Uranus", symbol=u"\u26E2", R=25362 * u.km)
Neptune = Body.from_parameters(
    Sun, k=6836529 * u.km ** 3 / u.s ** 2,
    name="Neptune", symbol=u"\u2646", R=24622 * u.km)
