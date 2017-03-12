# coding: utf-8
import pytest

from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from numpy import arange

from poliastro import bodies
from poliastro.twobody.orbit import Orbit

def test_body_has_k_given_in_constructor():
    k = 3.98e5 * u.km ** 3 / u.s ** 2
    earth = bodies.Body(k)
    assert earth.k == k


def test_body_from_parameters_raises_valueerror_if_k_units_not_correct():
    wrong_k = 4902.8 * u.kg
    _name = _symbol = ""
    _R = 0
    _parent = None
    _orbit = None
    with pytest.raises(u.UnitsError) as excinfo:
        moon = bodies.Body.from_parameters(wrong_k, _name, _symbol,
                                           _R, _parent, _orbit)
    assert ("UnitsError: Argument 'k' to function 'from_parameters' must be in units convertible to 'km3 / s2'."
            in excinfo.exconly())


def test_body_printing_has_name_and_symbol():
    name = "2 Pallas"
    symbol = u"\u26b4"
    k = 1.41e10 * u.m ** 3 / u.s ** 2
    pallas2 = bodies.Body(k, name, symbol)
    assert name in str(pallas2)
    assert symbol in str(pallas2)


def test_earth_has_k_given_in_literature():
    expected_k = 398600.44 * u.km ** 3 / u.s ** 2
    k = bodies.Earth.k
    assert_quantity_allclose(k.decompose([u.km, u.s]), expected_k)

def test_soi_calculation_working():
    _name = 'Planet'
    _k = 1 * u.km**3 / u.s**2
    _parent_name = 'Start'
    _parent_k = 32* _k
    _parent = bodies.Body.from_parameters( _parent_k, _parent_name, _parent_name,
                                         1 * u.km, None, None)
    _orbit = Orbit.circular(_parent, 3 * u.km)
    _planet = bodies.Body.from_parameters(_k, _name, _name,
                                          0.1 * u.km, _parent, _orbit)
    _planet.calculate_soi()
    _soi = _planet.soi
    expected_soi = 1 * u.km
    assert_quantity_allclose(_soi.to(u.km), expected_soi)
    
def test_generalized_Titius_Bode_Law():
    _bodies = ['Mercury', 'Venus', 'Earth', 'Mars', 'Ceres',
               'Jupiter', 'Saturn', 'Uranus', 'Neptune']
    au = 149597870.700 * u.km 
    expected_approx = au * 0.2075 * 1.7327 ** arange(1,10)
    obtained = []
    for _name in _bodies:
        obtained.append(bodies.body_dict[_name].orbit.a)
    assert_quantity_allclose(obtained, expected_approx,rtol = 0.2)
    
def test_main_bodies_soi_relative_magnitude():
    _bodies = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter',
               'Saturn', 'Uranus', 'Neptune', 'Moon']
    _expected = [46, 102, 145, 170, 687, 1025, 2040, 3525, 38]
    obtained = []
    for _name in _bodies:
        _body = bodies.body_dict[_name]
        obtained.append(_body.soi / _body.R)    
    assert_quantity_allclose(obtained, _expected,rtol = 0.2)
#bodies.