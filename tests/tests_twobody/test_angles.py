import astropy.time
import numpy as np
import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from numpy.testing import assert_allclose

from poliastro.bodies import Earth
from poliastro.core.elements import coe2mee, coe2rv, mee2coe, rv2coe
from poliastro.twobody import angles

# Data from Schlesinger & Udick, 1912
ELLIPTIC_ANGLES_DATA = [
    # ecc, M (deg), nu (deg)
    (0.0, 0.0, 0.0),
    (0.05, 10.0, 11.06),
    (0.06, 30.0, 33.67),
    (0.04, 120.0, 123.87),
    (0.14, 65.0, 80.50),
    (0.19, 21.0, 30.94),
    (0.35, 65.0, 105.71),
    (0.48, 180.0, 180.0),
    (0.75, 125.0, 167.57),
]


@pytest.fixture()
def classical():
    p = 11067.790  # u.km
    ecc = 0.83285  # u.one
    inc = np.deg2rad(87.87)  # u.rad
    raan = np.deg2rad(227.89)  # u.rad
    argp = np.deg2rad(53.38)  # u.rad
    nu = np.deg2rad(92.335)  # u.rad
    expected_res = (p, ecc, inc, raan, argp, nu)
    return expected_res


@pytest.fixture()
def circular():
    k = 3.9860047e14
    p = 24464560.0
    ecc = 0.0
    inc = 0.122138
    raan = 1.00681
    argp = 0.0
    nu = 0.048363
    expected_res = (p, ecc, inc, raan, argp, nu)
    return k, expected_res


@pytest.fixture()
def hyperbolic():
    k = 3.9860047e14
    p = 4.884856334147761e7
    ecc = 1.7311
    inc = 0.122138
    raan = 1.00681
    argp = 3.10686
    nu = 0.12741601769795755
    expected_res = (p, ecc, inc, raan, argp, nu)
    return k, expected_res


@pytest.fixture()
def equatorial():
    k = 3.9860047e14
    p = 1.13880762905224e7
    ecc = 0.7311
    inc = 0.0
    raan = 0.0
    argp = 3.10686
    nu = 0.44369564302687126
    expected_res = (p, ecc, inc, raan, argp, nu)
    return k, expected_res


@pytest.fixture()
def circular_equatorial():
    k = 3.9860047e14
    p = 1.13880762905224e7
    ecc = 0.0
    inc = 0.0
    raan = 0.0
    argp = 0.0
    nu = 0.44369564302687126
    expected_res = (p, ecc, inc, raan, argp, nu)
    return k, expected_res


def test_true_to_eccentric():
    # Data from NASA-TR-R-158
    data = [
        # ecc,E (deg), nu(deg)
        (0.0, 0.0, 0.0),
        (0.05, 10.52321, 11.05994),
        (0.10, 54.67466, 59.49810),
        (0.35, 142.27123, 153.32411),
        (0.61, 161.87359, 171.02189),
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
    for row in ELLIPTIC_ANGLES_DATA:
        ecc, M, expected_nu = row
        ecc = ecc * u.one
        M = M * u.deg
        expected_nu = expected_nu * u.deg

        nu = angles.E_to_nu(angles.M_to_E(M, ecc), ecc)

        assert_quantity_allclose(nu, expected_nu, rtol=1e-4)


def test_true_to_mean():
    for row in ELLIPTIC_ANGLES_DATA:
        ecc, expected_M, nu = row
        ecc = ecc * u.one
        expected_M = expected_M * u.deg
        nu = nu * u.deg

        M = angles.E_to_M(angles.nu_to_E(nu, ecc), ecc)

        assert_quantity_allclose(M, expected_M, rtol=1e-4)


def test_true_to_mean_hyperbolic():
    # Data from Curtis, H. (2013). *Orbital mechanics for engineering students*.
    # Example 3.5
    nu = 100 * u.deg
    ecc = 2.7696 * u.one
    expected_M = 11.279 * u.rad

    M = angles.F_to_M(angles.nu_to_F(nu, ecc), ecc)

    assert_quantity_allclose(M, expected_M, rtol=1e-4)


def test_mean_to_true_hyperbolic():
    # Data from Curtis, H. (2013). *Orbital mechanics for engineering students*.
    # Example 3.5
    M = 11.279 * u.rad
    ecc = 2.7696 * u.one
    expected_nu = 100 * u.deg

    nu = angles.F_to_nu(angles.M_to_F(M, ecc), ecc)

    assert_quantity_allclose(nu, expected_nu, rtol=1e-4)


def test_flight_path_angle():
    # Data from Curtis, example 2.5
    nu = 109.5 * u.deg
    ecc = 0.6 * u.one
    expected_gamma = 35.26 * u.deg

    gamma = angles.fp_angle(np.deg2rad(nu), ecc)

    assert_quantity_allclose(gamma, expected_gamma, rtol=1e-3)


@pytest.mark.parametrize(
    "expected_nu", np.linspace(-1 / 3.0, 1 / 3.0, num=100) * np.pi * u.rad
)
@pytest.mark.parametrize("ecc", [3200 * u.one, 1.5 * u.one])
def test_mean_to_true_hyperbolic_highecc(expected_nu, ecc):
    M = angles.F_to_M(angles.nu_to_F(expected_nu, ecc), ecc)
    nu = angles.F_to_nu(angles.M_to_F(M, ecc), ecc)
    assert_quantity_allclose(nu, expected_nu, rtol=1e-4)


@pytest.mark.parametrize("E", np.linspace(-1, 1, num=10) * np.pi * u.rad)
@pytest.mark.parametrize("ecc", np.linspace(0.1, 0.9, num=10) * u.one)
def test_eccentric_to_true_range(E, ecc):
    nu = angles.E_to_nu(E, ecc)
    E_result = angles.nu_to_E(nu, ecc)
    assert_quantity_allclose(E_result, E, rtol=1e-8)


def test_convert_between_coe_and_rv_is_transitive(classical):
    k = Earth.k.to(u.km ** 3 / u.s ** 2).value  # u.km**3 / u.s**2
    res = rv2coe(k, *coe2rv(k, *classical))
    assert_allclose(res, classical)


def test_convert_between_coe_and_mee_is_transitive(classical):
    res = mee2coe(*coe2mee(*classical))
    assert_allclose(res, classical)


def test_convert_coe_and_rv_circular(circular):
    k, expected_res = circular
    res = rv2coe(k, *coe2rv(k, *expected_res))
    assert_allclose(res, expected_res, atol=1e-8)


def test_convert_coe_and_rv_hyperbolic(hyperbolic):
    k, expected_res = hyperbolic
    res = rv2coe(k, *coe2rv(k, *expected_res))
    assert_allclose(res, expected_res, atol=1e-8)


def test_convert_coe_and_rv_equatorial(equatorial):
    k, expected_res = equatorial
    res = rv2coe(k, *coe2rv(k, *expected_res))
    assert_allclose(res, expected_res, atol=1e-8)


def test_convert_coe_and_rv_circular_equatorial(circular_equatorial):
    k, expected_res = circular_equatorial
    res = rv2coe(k, *coe2rv(k, *expected_res))
    assert_allclose(res, expected_res, atol=1e-8)


def test_raan_from_ltan_metopb():
    # MetOp-B LTAN: 21:31:45
    # LTAN from https://www.ospo.noaa.gov/Operations/METOP/status.html
    # METOP-B
    # 1 38771U 12049A   20049.95408566 -.00000014  00000-0  13607-4 0  9997
    # 2 38771  98.7092 110.9899 0001263  48.5458 295.8781 14.21485930385043

    ltan = (21 + ((31 + 45 / 60) / 60)) * u.hourangle
    epoch = astropy.time.Time(
        astropy.time.Time("2020-01-01 00:00").to_value("mjd") + 49.95408566 - 1,
        format="mjd",
    )
    expected_raan = 110.9899 * u.deg
    raan = angles.raan_from_ltan(epoch, ltan)
    assert_allclose(raan.wrap_at(360 * u.deg).to(u.deg), expected_raan, atol=0.3)


def test_raan_from_ltan_sentinel5p():
    # SENTINEL-5P LTAN: 13:30
    # LTAN from https://sentinels.copernicus.eu/web/sentinel/missions/sentinel-5p/geographical-coverage
    # 1 42969U 17064A   20049.78099017 -.00000032  00000-0  54746-5 0  9991
    # 2 42969  98.7249 350.5997 0001077  82.0109 278.1189 14.19549365121775

    ltan = (13 + (30 / 60)) * u.hourangle
    epoch = astropy.time.Time(
        astropy.time.Time("2020-01-01 00:00").to_value("mjd") + 49.78099017 - 1,
        format="mjd",
    )
    expected_raan = 350.5997 * u.deg
    raan = angles.raan_from_ltan(epoch, ltan)
    assert_allclose(raan.wrap_at(360 * u.deg).to(u.deg), expected_raan, atol=0.3)
