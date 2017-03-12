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

_data_path = os.path.join(os.path.dirname(__file__), "bodies_data.txt")
with open(_data_path, 'r', encoding='utf-8') as _raw_data :
    for _line in _raw_data:
        _line = _line.replace('\ufeff','')
        _line = _line.replace('\n','')
        if _line[0] == '#':
            continue
        _line = _line.replace(' ','')
        _body_data = _line.split(sep = ',')
        _name, _symbol, _k, _R = _body_data[:4]
        _k = float(_k) * u.km**3 / u.s**2
        _R = float(_R) * u.km
        
        if _body_data[4] == 'None':
            _parent = None
        else:
            _parent_name = _body_data[4]
            if _parent_name in body_dict:
                _parent = body_dict[_parent_name]
            else:
                _message = 'loading bodies warning: '
                _message += 'object not found during body loading:'
                _message += _parent_name +'. Parent reverted to None for ' + _name
                warnings.warn(_message)
                _parent = None
        if _body_data[5] == 'None':
            _orbit = None
        elif _parent == None:
            _orbit = None
            warnings.warn('Unable to create orbit without valid parent body')
        else:
            _a, _ecc, _inc, _L, _long_peri, _raan = _body_data[5:11]
            _a = float(_a) * au
            _ecc = float(_ecc) * u.one
            _inc = (float(_inc) * u.deg).to(u.rad)
            _L = (float(_L) * u.deg).to(u.rad)
            _long_peri = (float(_long_peri) * u.deg).to(u.rad)
            _raan = (float(_raan) * u.deg).to(u.rad)
            _argp = _long_peri - _raan
            _M = _L - _long_peri
            _nu = M_to_nu(_M, _ecc)
            _orbit = Orbit.from_classical(_parent, _a, _ecc, _inc, _raan, _argp, _nu)
        body_dict[_name] = Body.from_parameters(_k, _name, _symbol, _R, _parent, _orbit)
    

for _body_name in body_dict:
    _body = body_dict[_body_name]
    if _body.parent != None and _body.orbit != None:
        _body.calculate_soi()

Sun = body_dict['Sun']
Earth = body_dict['Earth']
Jupiter = body_dict['Jupiter']


