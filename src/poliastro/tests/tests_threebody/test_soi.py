import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from poliastro.bodies import (
    Earth,
    Jupiter,
    Mars,
    Mercury,
    Neptune,
    Saturn,
    Uranus,
    Venus,
)
from poliastro.threebody.soi import hill_radius, laplace_radius


@pytest.mark.parametrize(
    "body, expected_r_SOI",
    [
        (Mercury, 1.12e8 * u.m),
        (Venus, 6.16e8 * u.m),
        (Earth, 9.25e8 * u.m),
        (Mars, 5.77e8 * u.m),
        (Jupiter, 4.82e10 * u.m),
        (Saturn, 5.48e10 * u.m),
        (Uranus, 5.18e10 * u.m),
        (Neptune, 8.66e10 * u.m),
    ]
    # Data from Table A.2., Curtis "Orbital Mechanics for Engineering Students"
)
def test_laplace_radius(body, expected_r_SOI):
    r_SOI = laplace_radius(body)

    assert_quantity_allclose(r_SOI, expected_r_SOI, rtol=1e-1)


@pytest.mark.parametrize(
    "body, expected_r_SOI",
    [
        pytest.param(Mercury, 2.21e8 * u.m, marks=pytest.mark.xfail),  # Chebotarev
        (Mercury, 1.75e8 * u.m),  # Our result
        (Venus, 1.03e9 * u.m),
        (Earth, 1.49e9 * u.m),
        (Mars, 1.07e9 * u.m),
        (Jupiter, 5.28e10 * u.m),
        (Saturn, 6.50e10 * u.m),
        (Uranus, 7.01e10 * u.m),
        (Neptune, 1.16e11 * u.m),
    ],
    # Data from Chebotarev "Gravitational Spheres of the Major Planets, Moon and Sun",
    # notice the xfail for Mercury because we use the eccentricity
)
def test_hill_radius(body, expected_r_SOI):
    r_SOI = hill_radius(body)

    assert_quantity_allclose(r_SOI, expected_r_SOI, rtol=1e-1)
