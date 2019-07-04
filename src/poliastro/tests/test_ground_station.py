from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from poliastro.constants import J2000
from poliastro.ground_station import GroundStation


def test_cartesian_coordinates():
    expected_cords = [
        3764258.64785411 * u.m,
        3295359.33856106 * u.m,
        3942945.28570563 * u.m,
    ]

    el_cords = [38.43 * u.deg, 41.2 * u.deg, 0 * u.m]
    ellipsoid = [6378137 * u.m, 6356752.314245 * u.m]

    gs = GroundStation(*el_cords, J2000, 360.9856 * u.rad / u.day, *ellipsoid)
    c_cords = gs.cartesian_cords

    assert_quantity_allclose(c_cords, expected_cords)


def test_tangential_vectors():
    el_cords = [38.43 * u.deg, 41.2 * u.deg, 0 * u.m]
    ellipsoid = [6378137 * u.m, 6356752.314245 * u.m]

    gs = GroundStation(*el_cords, J2000, 360.9856 * u.rad / u.day, *ellipsoid)

    N = gs.N
    v1, v2 = gs.tangential_vecs

    assert abs(N.dot(v1)) <= 1e-7
    assert abs(N.dot(v2)) <= 1e-7


def test_visible():
    el_cords = [38.43 * u.deg, 41.2 * u.deg, 0 * u.m]
    ellipsoid = [6378137 * u.m, 6356752.314245 * u.m]

    gs = GroundStation(*el_cords, J2000, 360.9856 * u.rad / u.day, *ellipsoid)

    cords = gs.cartesian_cords

    p1 = [cords[i] + 10 * gs.N[i] * u.m for i in range(3)]
    p2 = [cords[i] - 10 * gs.N[i] * u.m for i in range(3)]

    # TODO: Fix assert failing with 'is'

    assert gs.is_visible(*p1) == True
    assert gs.is_visible(*p2) == False


def test_propagate():
    expected_cords = [
        3707017.86447794 * u.m,
        3359621.24213846 * u.m,
        3942945.2857615 * u.m,
    ]

    el_cords = [38.43 * u.deg, 41.2 * u.deg, 0 * u.m]
    ellipsoid = [6378137 * u.m, 6356752.314245 * u.m]

    gs = GroundStation(*el_cords, J2000, 360.9856 * u.deg / u.day, *ellipsoid)
    gs.propagate(J2000 + 1 * u.day)

    c_cords = gs.cartesian_cords

    assert_quantity_allclose(c_cords, expected_cords)
