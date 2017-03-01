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
    def __init__(self, k, name=None, symbol=None, R=0 * u.km):
        """Constructor.

        Parameters
        ----------
        k : Quantity
            Standard gravitational parameter.
        name : str
            Name of the body.
        symbol : str
            Symbol for the body.
        R : Quantity
            Radius of the body.

        """
        self.k = k
        self.name = name
        self.symbol = symbol
        self.R = R

    @classmethod
    @u.quantity_input(k=u.km ** 3 / u.s ** 2, R=u.km)
    def from_parameters(cls, k, name, symbol, R):
        return cls(k, name, symbol, R)

    def __str__(self):
        return u"{0} ({1})".format(self.name, self.symbol)

    def _repr_latex_(self):
        """Creates a LaTeX representation.

        Used by the IPython notebook.

        """
        return self.__str__()


Sun = Body.from_parameters(k=132712440018 * u.km ** 3 / u.s ** 2,
           name="Sun", symbol=u"\u2609", R=695700 * u.km)

Earth = Body.from_parameters(k=398600.44 * u.km ** 3 / u.s ** 2,
             name="Earth", symbol=u"\u2641", R=R_earth.to(u.km))
Earth.a = 149598023 * u.km #Semi-mayor axis
Earth.e = 0.0167086 * u.dimensionless_unscaled #Eccentricity

Moon = Body.from_parameters(k=4902.8 * u.km ** 3 / u.s ** 2,
             name="Moon", symbol=u"\u263E", R=1737.1 * u.km)
Moon.a = 384399 * u.km #Semi-mayor axis
Moon.e = 0.0549 * u.dimensionless_unscaled #Eccentricity

Mercury = Body.from_parameters(k=22032.09 * u.km ** 3 / u.s ** 2,
             name="Mercury", symbol=u"\u263F", R=2439.7 * u.km)
Mercury.a = 57909050 * u.km #Semi-mayor axis
Mercury.e = 0.205630 * u.dimensionless_unscaled #Eccentricity

Venus = Body.from_parameters(k=324858 * u.km ** 3 / u.s ** 2,
             name="Venus", symbol=u"\u2640", R=6051.8 * u.km)
Venus.a = 108208000 * u.km #Semi-mayor axis
Venus.e = 0.006772 * u.dimensionless_unscaled #Eccentricity

Mars = Body.from_parameters(k=42828.4 * u.km ** 3 / u.s ** 2,
             name="Mars", symbol=u"\u2642", R=3389.5 * u.km)
Mars.a = 227.9392 * u.Gm #Semi-mayor axis
Mars.e = 0.0934 * u.dimensionless_unscaled #Eccentricity

Jupiter = Body.from_parameters(k=126712764 * u.km ** 3 / u.s ** 2,
             name="Jupiter", symbol=u"\u2643", R=69911 * u.km)
Jupiter.a = 778.299 * u.Gm #Semi-mayor axis
Jupiter.e = 0.048498 * u.dimensionless_unscaled #Eccentricity

Saturn = Body.from_parameters(k=37940585 * u.km ** 3 / u.s ** 2,
             name="Saturn", symbol=u"\u2644", R=58232 * u.km)
Saturn.a = 1429.395 * u.Gm #Semi-mayor axis
Saturn.e = 0.05555 * u.dimensionless_unscaled #Eccentricity

Uranus = Body.from_parameters(k=5794548 * u.km ** 3 / u.s ** 2,
             name="Uranus", symbol=u"\u2645", R=25362 * u.km)
Uranus.a = 2875.04 * u.Gm #Semi-mayor axis
Uranus.e = 0.04638 * u.dimensionless_unscaled #Eccentricity

Neptune = Body.from_parameters(k=6836535 * u.km ** 3 / u.s ** 2,
             name="Neptune", symbol=u"\u2646", R=24622 * u.km)
Neptune.a = 4504.45 * u.Gm #Semi-mayor axis
Neptune.e = 0.009456 * u.dimensionless_unscaled #Eccentricity

Pluto = Body.from_parameters(k=977 * u.km ** 3 / u.s ** 2,
             name="Pluto", symbol=u"\u2647", R=1187 * u.km)
Pluto.a = 5906.38 * u.Gm #Semi-mayor axis
Pluto.e = 0.2488 * u.dimensionless_unscaled #Eccentricity

planet_list = [
        Mercury,
        Venus,
        Earth,
        Mars,
        Jupiter,
        Saturn,
        Uranus,
        Neptune,
        Pluto] #Kinda? Just so it doesn't feel alone

for planet in planet_list:
    planet.r_soi = planet.a * (planet.k / Sun.k) ** (2/5)   

Moon.r_soi = Moon.a * (Moon.k / Earth.k) ** (2/5)     


#Checking the numers:
if __name__ == '__main__':
    for planet in planet_list:
        print(planet.name.rjust(8), round(float(planet.r_soi / planet.R)))
    #Values should be close to the described at:
    #https://en.wikipedia.org/wiki/Sphere_of_influence_(astrodynamics)
