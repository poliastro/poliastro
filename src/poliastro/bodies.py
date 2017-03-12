# coding: utf-8
"""Bodies of the Solar System.

Contains some predefined bodies of the Solar System:

* Sun (☉)
* Earth (♁)

and a way to define new bodies (:py:class:`~Body` class).

"""

from astropy import units as u
import warnings
from poliastro.twobody import Orbit
from poliastro.twobody.angles import M_to_nu
import os.path

au = 149597870.700 * u.km
body_dict = {}

class Body(object):
    """Class to represent a body of the Solar System.

    """
    def __init__(self, k, name=None, symbol=None, R=0 * u.km,
                 parent=None, orbit = None):
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
        parent : Body
            Gravitational center (ej: Sun for the Earth)
        orbit: Orbit
            orbit object of the body

        """
        self.k = k
        self.name = name
        self.symbol = symbol
        self.R = R
        self.parent = parent
        self.orbit = orbit
        self.soi = None

    @classmethod
    @u.quantity_input(k=u.km ** 3 / u.s ** 2, R=u.km)
    def from_parameters(cls, k, name, symbol, R, parent, orbit):
        return cls(k, name, symbol, R, parent, orbit)

    def __str__(self):
        return u"{0} ({1})".format(self.name, self.symbol)

    def _repr_latex_(self):
        """Creates a LaTeX representation.

        Used by the IPython notebook.

        """
        return self.__str__()
    def calculate_soi(self):
        '''Calculate the radius of the sphere of influence.
        The body must have a parent body and an orbit defined'''
        if self.parent == None:
            warnings.warn('Unable to calculate SOI without valid parent body')
        elif self.orbit == None:
            warnings.warn('Unable to calculate SOI without valid orbit assigned')
        else:
            self.soi= self.orbit.a * (self.k / self.parent.k) ** (2/5)

data_path = os.path.join(os.path.dirname(__file__), "bodies_data.txt")
with open(data_path, 'r', encoding='utf-8') as raw_data :
    for line in raw_data:
        line = line.replace('\ufeff','')
        line = line.replace('\n','')
        if line[0] == '#':
            continue
        line = line.replace(' ','')
        body_data = line.split(sep = ',')
        name, symbol, k, R = body_data[:4]
        k = float(k) * u.km**3 / u.s**2
        R = float(R) * u.km
        
        if body_data[4] == 'None':
            parent = None
        else:
            parent_name = body_data[4]
            if parent_name in body_dict:
                parent = body_dict[parent_name]
            else:
                message = 'loading bodies warning: '
                message += 'object not found during body loading:'
                message += parent_name +'. Parent reverted to None for ' + name
                warnings.warn(message)
                parent = None
        if body_data[5] == 'None':
            orbit = None
        elif parent == None:
            orbit = None
            warnings.warn('Unable to create orbit without valid parent body')
        else:
            a, ecc, inc, L, long_peri, raan = body_data[5:11]
            a = float(a) * au
            ecc = float(ecc) * u.one
            inc = (float(inc) * u.deg).to(u.rad)
            L = (float(L) * u.deg).to(u.rad)
            long_peri = (float(long_peri) * u.deg).to(u.rad)
            raan = (float(raan) * u.deg).to(u.rad)
            argp = long_peri - raan
            M = L - long_peri
            nu = M_to_nu(M, ecc)
            orbit = Orbit.from_classical(parent, a, ecc, inc, raan, argp, nu)
        body_dict[name] = Body.from_parameters(k, name, symbol, R, parent, orbit)
    

for body_name in body_dict:
    body = body_dict[body_name]
    if body.parent != None and body.orbit != None:
        body.calculate_soi()

Sun = body_dict['Sun']
Earth = body_dict['Earth']
Jupiter = body_dict['Jupiter']

#Checking the numers:
if __name__ == '__main__':
    for body_name in body_dict:
        body = body_dict[body_name]
        if body.soi != None:
            print(body.name.rjust(8), round(float(body.soi / body.R)))
    #Values should be close to the described at:
    #https://en.wikipedia.org/wiki/Sphere_of_influence_(astrodynamics)
