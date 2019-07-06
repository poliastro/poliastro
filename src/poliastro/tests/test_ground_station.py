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


def test_tangential_vectors():
    el_cords = [38.43 * u.deg, 41.2 * u.deg, 0 * u.m]
    ellipsoid = [6378137 * u.m, 6356752.314245 * u.m]

    gs = GroundStation(*el_cords, *ellipsoid)

    N = gs.N
    v1, v2 = gs.tangential_vecs

    assert abs(N.dot(v1)) <= 1e-7
    assert abs(N.dot(v2)) <= 1e-7


def test_visible():
    el_cords = [38.43 * u.deg, 41.2 * u.deg, 0 * u.m]
    ellipsoid = [6378137 * u.m, 6356752.314245 * u.m]

    gs = GroundStation(*el_cords, *ellipsoid)

    cords = gs.cartesian_cords

    p1 = [cords[i] + 10 * gs.N[i] * u.m for i in range(3)]
    p2 = [cords[i] - 10 * gs.N[i] * u.m for i in range(3)]

    assert gs.is_visible(*p1)
    assert not gs.is_visible(*p2)


def test_cartesian_conversion_approximate():
    el_cords = [0.670680 * u.rad, 0.7190227 * u.rad, 0 * u.m]

    c_cords = [3764258.64785411 * u.m, 3295359.33856106 * u.m, 3942945.28570563 * u.m]
    ellipsoid = [6378137 * u.m, 6356752.314245 * u.m]

    gs = GroundStation(*el_cords, *ellipsoid)

    cords = gs.cartesian_to_ellipsoidal(*c_cords)

    # TODO: Find ways to improve error margin
    assert_quantity_allclose(el_cords[0], cords[0], 1e-4)
    assert_quantity_allclose(el_cords[1], cords[1], 1e-4)
    assert_quantity_allclose(el_cords[2], cords[2], 1)
