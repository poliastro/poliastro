# coding: utf-8
import pytest
from astropy import units as u

from poliastro import bodies


def test_body_has_k_given_in_constructor():
    k = 3.98e5 * u.km ** 3 / u.s ** 2
    earth = bodies.Body(k)
    assert earth.k == k


def test_body_constructor_raises_valueerror_if_k_units_not_correct():
    wrong_k = 4902.8 * u.kg
    with pytest.raises(u.UnitsError) as excinfo:
        moon = bodies.Body(wrong_k)
    assert ("UnitsError: Units must be consistent"
            in excinfo.exconly())


def test_body_printing_has_name_and_symbol():
    name = "2 Pallas"
    symbol = u"\u26b4"
    k = 1.41e10 * u.m ** 3 / u.s ** 2
    pallas2 = bodies.Body(k, name, symbol)
    # NOTE: str(pallas2) fails in Python 2, not willing to provide a fix; see
    # http://docs.python.org/3/howto/pyporting.html#str-unicode
    assert name in pallas2.__str__()
    assert symbol in pallas2.__str__()
