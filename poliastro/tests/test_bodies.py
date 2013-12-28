import pytest
from astropy import units as u

from poliastro import bodies


def test_body_has_k_given_in_constructor():
    k = 3.98e5 * u.km ** 3 / u.s ** 2
    earth = bodies.Body(k)
    assert earth.k == k


def test_body_constructor_raises_valueerror_if_k_is_not_quantity():
    k = 4902.8
    with pytest.raises(ValueError) as excinfo:
        moon = bodies.Body(k)
    assert ("ValueError: k must have units (use astropy.units)"
            in excinfo.exconly())


def test_body_constructor_raises_valueerror_if_k_units_not_correct():
    k = 4902.8 * u.kg
    with pytest.raises(ValueError) as excinfo:
        moon = bodies.Body(k)
    assert ("k units not consistent (expected u.m ** 3 / u.s ** 2)"
            in excinfo.exconly())


def test_body_printing_has_name_and_symbol():
    name = "2 Pallas"
    symbol = u"\u26b4"
    k = 1.41e10 * u.m ** 3 / u.s ** 2
    pallas2 = bodies.Body(k, name, symbol)
    assert name in str(pallas2)
    assert symbol in str(pallas2)
