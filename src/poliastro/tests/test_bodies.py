import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from poliastro import bodies


def test_body_has_k_given_in_constructor():
    k = 3.98e5 * u.km ** 3 / u.s ** 2
    earth = bodies.Body(None, k, "")
    assert earth.k == k


def test_body_from_parameters_raises_valueerror_if_k_units_not_correct():
    wrong_k = 4902.8 * u.kg
    _name = _symbol = ""
    _R = 0
    with pytest.raises(u.UnitsError) as excinfo:
        bodies.Body.from_parameters(None, wrong_k, _name, _symbol, _R)
    assert (
        "UnitsError: Argument 'k' to function 'from_parameters' must be in units convertible to 'km3 / s2'."
        in excinfo.exconly()
    )


def test_body_printing_has_name_and_symbol():
    name = "2 Pallas"
    symbol = u"\u26b4"
    k = 1.41e10 * u.m ** 3 / u.s ** 2
    pallas2 = bodies.Body(None, k, name, symbol)
    assert name in str(pallas2)
    assert symbol in str(pallas2)


def test_earth_has_k_given_in_literature():
    expected_k = 3.986004418e14 * u.m ** 3 / u.s ** 2
    k = bodies.Earth.k
    assert_quantity_allclose(k.decompose([u.km, u.s]), expected_k)


def test_body_kwargs():
    name = "2 Pallas"
    symbol = u"\u26b4"
    k = 1.41e10 * u.m ** 3 / u.s ** 2
    pallas2 = bodies.Body(None, k, name, symbol)
    assert pallas2.kwargs == {}
    pallas2 = bodies.Body(None, k, name, symbol, custom="data")
    assert "custom" in pallas2.kwargs


def test_from_relative():
    TRAPPIST1 = bodies.Body.from_relative(
        reference=bodies.Sun,
        parent=None,
        k=0.08,  # Relative to the Sun
        name="TRAPPIST",
        symbol=None,
        R=0.114,
    )  # Relative to the Sun

    # check values properly calculated
    VALUECHECK = bodies.Body.from_relative(
        reference=bodies.Earth,
        parent=TRAPPIST1,
        k=1,
        name="VALUECHECK",
        symbol=None,
        R=1,
    )
    assert bodies.Earth.k == VALUECHECK.k
    assert bodies.Earth.R == VALUECHECK.R
