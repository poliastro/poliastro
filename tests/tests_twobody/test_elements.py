from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from poliastro.twobody.elements import circular_velocity


def test_simple_circular_velocity():
    k = 398600 * u.km ** 3 / u.s ** 2
    a = 7000 * u.km

    expected_V = 7.5460491 * u.km / u.s

    V = circular_velocity(k, a)

    assert_quantity_allclose(V, expected_V)
