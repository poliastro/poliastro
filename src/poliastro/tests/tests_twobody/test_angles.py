import numpy as np

from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

import pytest
from poliastro.twobody import angles


def test_true_to_eccentric():
    # Data from NASA-TR-R-158
    data = [
        # ecc,E (deg), nu(deg)
        (0.0, 0.0, 0.0),
        (0.05, 10.52321, 11.05994),
        (0.10, 54.67466, 59.49810),
        (0.35, 142.27123, 153.32411),
        (0.61, 161.87359, 171.02189)
    ]
    for row in data:
        ecc, expected_E, nu = row
        ecc = ecc * u.one
        expected_E = expected_E * u.deg
        nu = nu * u.deg

        E = angles.nu_to_E(nu, ecc)

        assert_quantity_allclose(E, expected_E, rtol=1e-6)


def test_true_to_eccentric_hyperbolic():
    # Data from Curtis, H. (2013). *Orbital mechanics for engineering students*.
    # Example 3.5
    nu = 100 * u.deg
    ecc = 2.7696 * u.one
    expected_F = 2.2927 * u.rad

    F = angles.nu_to_F(nu, ecc)

    assert_quantity_allclose(F, expected_F, rtol=1e-4)


def test_mean_to_true():
    # Data from Schlesinger & Udick, 1912
    data = [
        # ecc, M (deg), nu (deg)
        (0.0, 0.0, 0.0),
        (0.05, 10.0, 11.06),
        (0.06, 30.0, 33.67),
        (0.04, 120.0, 123.87),
        (0.14, 65.0, 80.50),
        (0.19, 21.0, 30.94),
        (0.35, 65.0, 105.71),
        (0.48, 180.0, 180.0),
        (0.75, 125.0, 167.57)
    ]
    for row in data:
        ecc, M, expected_nu = row
        ecc = ecc * u.one
        M = M * u.deg
        expected_nu = expected_nu * u.deg

        nu = angles.M_to_nu(M, ecc)

        assert_quantity_allclose(nu, expected_nu, rtol=1e-4)


def test_true_to_mean():
    # Data from Schlesinger & Udick, 1912
    data = [
        # ecc, M (deg), nu (deg)
        (0.0, 0.0, 0.0),
        (0.05, 10.0, 11.06),
        (0.06, 30.0, 33.67),
        (0.04, 120.0, 123.87),
        (0.14, 65.0, 80.50),
        (0.19, 21.0, 30.94),
        (0.35, 65.0, 105.71),
        (0.48, 180.0, 180.0),
        (0.75, 125.0, 167.57)
    ]
    for row in data:
        ecc, expected_M, nu = row
        ecc = ecc * u.one
        expected_M = expected_M * u.deg
        nu = nu * u.deg

        M = angles.nu_to_M(nu, ecc)

        assert_quantity_allclose(M, expected_M, rtol=1e-4)


def test_true_to_mean_hyperbolic():
    # Data from Curtis, H. (2013). *Orbital mechanics for engineering students*.
    # Example 3.5
    nu = 100 * u.deg
    ecc = 2.7696 * u.one
    expected_M = 11.279 * u.rad

    M = angles.nu_to_M(nu, ecc)

    assert_quantity_allclose(M, expected_M, rtol=1e-4)


def test_mean_to_true_hyperbolic():
    # Data from Curtis, H. (2013). *Orbital mechanics for engineering students*.
    # Example 3.5
    M = 11.279 * u.rad
    ecc = 2.7696 * u.one
    expected_nu = 100 * u.deg

    nu = angles.M_to_nu(M, ecc)

    assert_quantity_allclose(nu, expected_nu, rtol=1e-4)


def test_flight_path_angle():
    # Data from Curtis, example 2.5
    nu = 109.5 * u.deg
    ecc = 0.6 * u.one
    expected_gamma = 35.26 * u.deg

    gamma = angles.fp_angle(np.deg2rad(nu), ecc)

    assert_quantity_allclose(gamma, expected_gamma, rtol=1e-3)


@pytest.mark.parametrize("expected_nu", np.linspace(-1 / 3.0, 1 / 3.0, num=100) * np.pi * u.rad)
@pytest.mark.parametrize("ecc", [3200 * u.one, 1.5 * u.one])
def test_mean_to_true_hyperbolic_highecc(expected_nu, ecc):
    M = angles.nu_to_M(expected_nu, ecc)
    print(M, ecc, M.value, ecc.value)
    nu = angles.M_to_nu(M, ecc)
    assert_quantity_allclose(nu, expected_nu, rtol=1e-4)
