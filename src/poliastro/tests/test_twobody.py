# coding: utf-8
import numpy as np
from numpy.testing import assert_array_almost_equal

from astropy import units as u

from poliastro.bodies import Earth

from poliastro.twobody.rv import rv2coe
from poliastro.twobody.classical import coe2rv, coe2mee
from poliastro.twobody.equinoctial import mee2coe


def test_convert_between_coe_and_rv_is_transitive():
    k = Earth.k.to(u.km**3 / u.s**2).value  # u.km**3 / u.s**2
    p = 11067.790  # u.km
    ecc = 0.83285  # u.one
    inc = np.deg2rad(87.87)  # u.rad
    raan = np.deg2rad(227.89)  # u.rad
    argp = np.deg2rad(53.38)  # u.rad
    nu = np.deg2rad(92.335)  # u.rad

    expected_res = (p, ecc, inc, raan, argp, nu)

    res = rv2coe(k, *coe2rv(k, *expected_res))

    assert_array_almost_equal(res, expected_res)


def test_convert_between_coe_and_mee_is_transitive():
    p = 11067.790  # u.km
    ecc = 0.83285  # u.one
    inc = np.deg2rad(87.87)  # u.rad
    raan = np.deg2rad(227.89)  # u.rad
    argp = np.deg2rad(53.38)  # u.rad
    nu = np.deg2rad(92.335)  # u.rad

    expected_res = (p, ecc, inc, raan, argp, nu)

    res = mee2coe(*coe2mee(*expected_res))

    assert_array_almost_equal(res, expected_res)
