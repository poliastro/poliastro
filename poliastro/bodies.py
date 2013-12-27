"""Bodies of the Solar System.

Contains some predefined bodies of the Solar System:

* Sun (☉)
* Earth (♁)

and a way to define new bodies (`Body` class).

"""

from astropy.constants import G, M_sun, M_earth
from astropy import units as u


class Body(object):
    """Class to represent a body of the Solar System.

    """
    def __init__(self, k, name=None, symbol=None):
        """Constructor.

        Parameters
        ----------
        k : Quantity
            Standard gravitational parameter

        """
        try:
            if k.si.unit != u.m ** 3 / u.s ** 2:
                raise ValueError("k units not consistent "
                                 "(expected u.m ** 3 / u.s ** 2)")
        except AttributeError:
            raise ValueError("k must have units (use astropy.units)")
        self.k = k
        self.name = name
        self.symbol = symbol

    def __str__(self):
        return "{} ({})".format(self.name, self.symbol)


Sun = Body(k=G * M_sun, name="Sun", symbol=u"☉")
Earth = Body(k = G * M_earth, name="Earth", symbol=u"♁")
