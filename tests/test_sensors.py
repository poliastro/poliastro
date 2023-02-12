from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
import pytest

from poliastro.bodies import Earth
from poliastro.sensors import (
    ground_range_diff_at_azimuth,
    min_and_max_ground_range,
)


# Example taken from "Fundamentals of Astrodynamics and Applications", 4th ed (2013)" by David A. Vallado, pages 859-860
@pytest.mark.parametrize(
    "altitude, fov, boresight, expected_lat_lon_max, expected_lat_lon_min",
    [
        (
            800 * u.km,
            (25 * u.deg).to(u.rad),
            (40 * u.deg).to(u.rad),
            10.73517 * u.deg,
            3.80977 * u.deg,
        ),
        (
            800 * u.km,
            (25 * u.deg).to(u.rad),
            (0 * u.deg).to(u.rad),
            1.5984 * u.deg,
            -1.5984 * u.deg,
        ),
    ],
)
def test_max_and_min_ground_range(
    altitude, fov, boresight, expected_lat_lon_max, expected_lat_lon_min
):

    R = Earth.R.to(u.km)
    lat_lon_min, lat_lon_max = min_and_max_ground_range(
        altitude, fov, boresight, R
    )
    assert_quantity_allclose(lat_lon_max, expected_lat_lon_max, rtol=1e-4)
    assert_quantity_allclose(lat_lon_min, expected_lat_lon_min, rtol=1e-4)


# Example taken from "Fundamentals of Astrodynamics and Applications", 4th ed (2013)" by David A. Vallado, pages 859-860
@pytest.mark.parametrize(
    "altitude, fov, boresight, azimuth, nadir_lat, nadir_lon, expected_ground_range_diff, expected_target_lat, expected_target_lon",
    [
        (
            800 * u.km,
            (25 * u.deg).to(u.rad),
            (40 * u.deg).to(u.rad),
            (140 * u.deg).to(u.rad),
            (50 * u.deg).to(u.rad),
            (40 * u.deg).to(u.rad),
            (6.9254 / 2) * u.deg,
            44.9926 * u.deg,
            45.7577 * u.deg,
        ),
    ],
)
def test_ground_range_diff_at_azimuth(
    altitude,
    fov,
    boresight,
    azimuth,
    nadir_lat,
    nadir_lon,
    expected_ground_range_diff,
    expected_target_lat,
    expected_target_lon,
):

    R = Earth.R.to(u.km)
    ground_range_diff, target_lat, target_lon = ground_range_diff_at_azimuth(
        altitude, fov, boresight, azimuth, nadir_lat, nadir_lon, R
    )
    assert_quantity_allclose(
        ground_range_diff, expected_ground_range_diff, rtol=1e-5
    )
    assert_quantity_allclose(target_lat, expected_target_lat, rtol=1e-6)
    assert_quantity_allclose(target_lon, expected_target_lon, rtol=1e-6)


# Example taken from "Fundamentals of Astrodynamics and Applications", 4th ed (2013)" by David A. Vallado, pages 859-860
@pytest.mark.parametrize(
    "altitude, fov, boresight, azimuth, nadir_lat, nadir_lon",
    [
        (
            800 * u.km,
            (25 * u.deg).to(u.rad),
            (40 * u.deg).to(u.rad),
            (190 * u.deg).to(u.rad),
            (50 * u.deg).to(u.rad),
            (40 * u.deg).to(u.rad),
        ),
    ],
)
def test_exception_ground_range_diff_at_azimuth(
    altitude,
    fov,
    boresight,
    azimuth,
    nadir_lat,
    nadir_lon,
):

    R = Earth.R.to(u.km)
    with pytest.raises(ValueError) as excinfo:
        ground_range_diff_at_azimuth(
            altitude, fov, boresight, azimuth, nadir_lat, nadir_lon, R
        )
    assert "beta must be between 0º and 180º" in excinfo.exconly()
