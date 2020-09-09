import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from poliastro.bodies import Earth
from poliastro.earth.sensors import (
    min_and_max_ground_range,
    max_and_min_ground_range_with_specific_azimuth,
)


# Example taken from "Fundamentals of Astrodynamics and Applications", 4th ed (2013)" by David A. Vallado, pages 859-860
@pytest.mark.parametrize(
    "h, n_fov, n_center,expected_Λ_max,expected_Λ_min ",
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
def test_max_and_min_ground_range(h, n_fov, n_center, expected_Λ_max, expected_Λ_min):

    R = Earth.R.to(u.km)
    Λ_min, Λ_max = min_and_max_ground_range(h, n_fov, n_center, R)
    assert_quantity_allclose(Λ_max, expected_Λ_max)
    assert_quantity_allclose(Λ_min, expected_Λ_min)


# Example taken from "Fundamentals of Astrodynamics and Applications", 4th ed (2013)" by David A. Vallado, pages 859-860
@pytest.mark.parametrize(
    "h, n_fov, n_center, β, φ_nadir, λ_nadir, expected_delta_Λ, expected_φ_tgt, expected_λ_tgt",
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
def test_max_and_min_ground_range_with_specific_azimuth(
    h,
    n_fov,
    n_center,
    β,
    φ_nadir,
    λ_nadir,
    expected_delta_Λ,
    expected_φ_tgt,
    expected_λ_tgt,
):

    R = Earth.R.to(u.km)
    delta_Λ, φ_tgt, λ_tgt = max_and_min_ground_range_with_specific_azimuth(
        h, n_center, n_fov, β, φ_nadir, λ_nadir, R
    )
    assert_quantity_allclose(delta_Λ, expected_delta_Λ)
    assert_quantity_allclose(φ_tgt, expected_φ_tgt)
    assert_quantity_allclose(λ_tgt, expected_λ_tgt)


# Example taken from "Fundamentals of Astrodynamics and Applications", 4th ed (2013)" by David A. Vallado, pages 859-860
@pytest.mark.parametrize(
    "h, n_fov, n_center, β, φ_nadir, λ_nadir, expected_delta_Λ, expected_φ_tgt, expected_λ_tgt",
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
def test_exception_max_and_min_ground_range_with_specific_azimuth(
    h,
    n_fov,
    n_center,
    β,
    φ_nadir,
    λ_nadir,
    expected_delta_Λ,
    expected_φ_tgt,
    expected_λ_tgt,
):

    R = Earth.R.to(u.km)
    with pytest.raises(ValueError) as excinfo:
        max_and_min_ground_range_with_specific_azimuth(
            h, n_center, n_fov, β, φ_nadir, λ_nadir, R
        )
    assert "β must be between 0º and 180º" in excinfo.exconly()
