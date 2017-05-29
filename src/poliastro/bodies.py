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
Earth = Body.from_parameters(
    Sun, k=398600 * u.km ** 3 / u.s ** 2,
    name="Earth", symbol=u"\u2641", R=R_earth.to(u.km))
