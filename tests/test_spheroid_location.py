from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from poliastro.bodies import Earth
from poliastro.spheroid_location import SpheroidLocation


def test_cartesian_coordinates():
    expected_cords = [
        3764258.64785411 * u.m,
        3295359.33856106 * u.m,
        3942945.28570563 * u.m,
    ]

    el_cords = (38.43 * u.deg, 41.2 * u.deg, 0 * u.m)

    p = SpheroidLocation(*el_cords, Earth)
    c_cords = p.cartesian_cords

    assert_quantity_allclose(c_cords, expected_cords)


def test_tangential_vectors():
    el_cords = (38.43 * u.deg, 41.2 * u.deg, 0 * u.m)

    p = SpheroidLocation(*el_cords, Earth)

    N = p.N
    v1, v2 = p.tangential_vecs

    assert abs(N.dot(v1)) <= 1e-7
    assert abs(N.dot(v2)) <= 1e-7


def test_visible():
    el_cords = (38.43 * u.deg, 41.2 * u.deg, 0 * u.m)

    p = SpheroidLocation(*el_cords, Earth)

    cords = p.cartesian_cords

    p1 = [cords[i] + 10 * p.N[i] * u.m for i in range(3)]
    p2 = [cords[i] - 10 * p.N[i] * u.m for i in range(3)]

    assert p.is_visible(*p1)
    assert not p.is_visible(*p2)


def test_f():
    expected_f = 0.0033528131

    el_cords = (38.43 * u.deg, 41.2 * u.deg, 0 * u.m)

    p = SpheroidLocation(*el_cords, Earth)

    f = p.f

    assert_quantity_allclose(f, expected_f)


def test_radius_of_curvature():
    expected_roc = 6363141.421601379 * u.m

    el_cords = (38.43 * u.deg, 41.2 * u.deg, 0 * u.m)

    p = SpheroidLocation(*el_cords, Earth)

    roc = p.radius_of_curvature

    assert_quantity_allclose(roc, expected_roc)


def test_distance():
    expected_distance = 6369864.745418392 * u.m
    el_cords = (38.43 * u.deg, 41.2 * u.deg, 0 * u.m)
    point_cords = (10.5 * u.m, 35.5 * u.m, 45.5 * u.m)

    p = SpheroidLocation(*el_cords, Earth)

    distance = p.distance(*point_cords)

    assert_quantity_allclose(distance, expected_distance)


def test_cartesian_conversion_approximate():
    el_cords = (0.670680 * u.rad, 0.7190227 * u.rad, 0 * u.m)

    c_cords = [3764258.64785411 * u.m, 3295359.33856106 * u.m, 3942945.28570563 * u.m]

    p = SpheroidLocation(*el_cords, Earth)

    cords = p.cartesian_to_ellipsoidal(*c_cords)

    # TODO: Find ways to improve error margin
    assert_quantity_allclose(el_cords[0], cords[0], 1e-4)
    assert_quantity_allclose(el_cords[1], cords[1], 1e-4)
    assert_quantity_allclose(el_cords[2], cords[2], 1)
