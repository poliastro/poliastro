# coding: utf-8
"""Bodies of the Solar System.

Contains some predefined bodies of the Solar System:

* Sun (☉)
* Earth (♁)

and a way to define new bodies (`Body` class).

"""

from astropy.constants import G, M_sun, M_earth, R_earth
from astropy import units as u

from poliastro.util import check_units


class Body(object):
    """Class to represent a body of the Solar System.

    """
    def __init__(self, k, name=None, symbol=None, R=None):
        """Constructor.

        Parameters
        ----------
        k : Quantity
            Standard gravitational parameter

        """
        if not check_units((k,), (u.m ** 3 / u.s ** 2,)):
            raise u.UnitsError("Units must be consistent")

        self.k = k
        self.name = name
        self.symbol = symbol
        self.R = R

    def __str__(self):
        return u"{} ({})".format(self.name, self.symbol)


Sun = Body(k=(G * M_sun).to(u.km ** 3 / u.s ** 2),
           name="Sun", symbol=u"\u2609")
Earth = Body(k=(G * M_earth).to(u.km ** 3 / u.s ** 2),
             name="Earth", symbol=u"\u2641", R=R_earth)
