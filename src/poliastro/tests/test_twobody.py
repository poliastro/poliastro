import pytest
import numpy as np
from numpy.testing import assert_allclose

from astropy import units as u

from poliastro.bodies import Earth
from poliastro.core.elements import coe2rv, coe2mee, rv2coe
from poliastro.twobody.equinoctial import mee2coe


@pytest.fixture()
def expected_res():
    p = 11067.790  # u.km
    ecc = 0.83285  # u.one
    inc = np.deg2rad(87.87)  # u.rad
    raan = np.deg2rad(227.89)  # u.rad
    argp = np.deg2rad(53.38)  # u.rad
    nu = np.deg2rad(92.335)  # u.rad
    expected_res = (p, ecc, inc, raan, argp, nu)
    return expected_res


def test_convert_between_coe_and_rv_is_transitive(expected_res):
    k = Earth.k.to(u.km**3 / u.s**2).value  # u.km**3 / u.s**2
    res = rv2coe(k, *coe2rv(k, *expected_res))
    assert_allclose(res, expected_res)


def test_convert_between_coe_and_mee_is_transitive(expected_res):
    res = mee2coe(*coe2mee(*expected_res))
    assert_allclose(res, expected_res)
