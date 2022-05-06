import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from poliastro.bodies import Earth
from poliastro.earth.sensors import (
    ground_range_diff_at_azimuth,
    min_and_max_ground_range,
)


# Example taken from "Fundamentals of Astrodynamics and Applications", 4th ed (2013)" by David A. Vallado, pages 859-860
@pytest.mark.parametrize(
    "h, eta_fov, eta_center, expected_lambda_max, expected_lambda_min ",
    [
        (
            800 * u.km,
            (25 * u.deg).to(u.rad),
            (40 * u.deg).to(u.rad),
            0.18736414 * u.rad,
            0.06649331 * u.rad,
        ),
        (
            800 * u.km,
            (25 * u.deg).to(u.rad),
            (0 * u.deg).to(u.rad),
            0.0278967 * u.rad,
            -0.0278967 * u.rad,
        ),
    ],
)
def test_max_and_min_ground_range(
    h, eta_fov, eta_center, expected_lambda_max, expected_lambda_min
):

    R = Earth.R.to(u.km)
    lambda_min, lambda_max = min_and_max_ground_range(
        h, eta_fov, eta_center, R
    )
    assert_quantity_allclose(lambda_max, expected_lambda_max)
    assert_quantity_allclose(lambda_min, expected_lambda_min)


# Example taken from "Fundamentals of Astrodynamics and Applications", 4th ed (2013)" by David A. Vallado, pages 859-860
@pytest.mark.parametrize(
    "h, eta_fov, eta_center, beta, phi_nadir, lambda_nadir, expected_delta_lambda, expected_phi_tgt, expected_lambda_tgt",
    [
        (
            800 * u.km,
            (25 * u.deg).to(u.rad),
            (40 * u.deg).to(u.rad),
            (140 * u.deg).to(u.rad),
            (50 * u.deg).to(u.rad),
            (40 * u.deg).to(u.rad),
            0.06197359 * u.rad,
            0.82639242 * u.rad,
            0.75442901 * u.rad,
        ),
    ],
)
def test_ground_range_diff_at_azimuth(
    h,
    eta_fov,
    eta_center,
    beta,
    phi_nadir,
    lambda_nadir,
    expected_delta_lambda,
    expected_phi_tgt,
    expected_lambda_tgt,
):

    R = Earth.R.to(u.km)
    delta_lambda, phi_tgt, lambda_tgt = ground_range_diff_at_azimuth(
        h, eta_center, eta_fov, beta, phi_nadir, lambda_nadir, R
    )
    assert_quantity_allclose(delta_lambda, expected_delta_lambda)
    assert_quantity_allclose(phi_tgt, expected_phi_tgt)
    assert_quantity_allclose(lambda_tgt, expected_lambda_tgt)


# Example taken from "Fundamentals of Astrodynamics and Applications", 4th ed (2013)" by David A. Vallado, pages 859-860
@pytest.mark.parametrize(
    "h, eta_fov, eta_center, beta, phi_nadir, lambda_nadir, expected_delta_lambda, expected_phi_tgt, expected_lambda_tgt",
    [
        (
            800 * u.km,
            (25 * u.deg).to(u.rad),
            (40 * u.deg).to(u.rad),
            (190 * u.deg).to(u.rad),
            (50 * u.deg).to(u.rad),
            (40 * u.deg).to(u.rad),
            0.06197359 * u.rad,
            0.82639242 * u.rad,
            0.75442901 * u.rad,
        ),
    ],
)
def test_exception_ground_range_diff_at_azimuth(
    h,
    eta_fov,
    eta_center,
    beta,
    phi_nadir,
    lambda_nadir,
    expected_delta_lambda,
    expected_phi_tgt,
    expected_lambda_tgt,
):

    R = Earth.R.to(u.km)
    with pytest.raises(ValueError) as excinfo:
        ground_range_diff_at_azimuth(
            h, eta_center, eta_fov, beta, phi_nadir, lambda_nadir, R
        )
    assert "beta must be between 0º and 180º" in excinfo.exconly()
