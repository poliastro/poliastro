from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from poliastro.ground_station import GroundStation


def test_cartesian_coordinates():
    expected_cords = [
        3764258.64785411 * u.m,
        3295359.33856106 * u.m,
        3942945.28570563 * u.m,
    ]

    el_cords = [38.43 * u.deg, 41.2 * u.deg, 0 * u.m]
    ellipsoid = [6378137 * u.m, 6356752.314245 * u.m]

    gs = GroundStation(*el_cords, *ellipsoid)
    c_cords = gs.cartesian_cords

    assert_quantity_allclose(c_cords, expected_cords)
