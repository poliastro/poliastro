# coding: utf-8
import pytest

from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from poliastro import bodies


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
