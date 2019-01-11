import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from poliastro.bodies import (
    Body,
    Earth,
    Jupiter,
    Mars,
    Mercury,
    Moon,
    Neptune,
    Pluto,
    Saturn,
    Sun,
    Uranus,
    Venus,
)
from poliastro.threebody.soi import laplace_radius, hill_radius


@pytest.mark.parametrize(
    "body, expected_r_SOI",
    [
        (Sun, None),
        (Mercury, 1.12e8),
        (Venus, 6.16e8),
        (Earth, 9.25e8),
        (Mars, 5.77e8),
        (Jupiter, 4.82e10),
        (Saturn, 5.48e10),
        (Uranus, 5.18e10),
        (Neptune, 8.66e10),
    ]
    # Data from Table A.2., Curtis "Orbital Mechanics for Engineering Students"
)
def test_laplace_radius(body, expected_r_SOI):
    if expected_r_SOI is not None:
        expected_r_SOI = expected_r_SOI * u.m

    r_SOI = laplace_radius(body)

    assert_quantity_allclose(r_SOI, expected_r_SOI, rtol=1e-1)


@pytest.mark.parametrize("missing_body", [Moon, Pluto])
def test_compute_missing_body_laplace_radius_raises_error(missing_body):
    with pytest.raises(RuntimeError) as excinfo:
        laplace_radius(missing_body)
    assert (
        "To compute the semimajor axis or eccentricity for Moon and Pluto use the JPL ephemeris"
        in excinfo.exconly()
    )


def test_laplace_radius_given_a():
    parent = Body(None, 1 * u.km ** 3 / u.s ** 2, "Parent")
    body = Body(parent, 1 * u.km ** 3 / u.s ** 2, "Body")
    r_SOI = laplace_radius(body, 1 * u.km)

    assert r_SOI == 1 * u.km


@pytest.mark.parametrize(
    "body, expected_r_SOI",
    [
        (Sun, None),
        (Mercury, 1.8e8),
        (Venus, 1.01e9),
        (Earth, 1.46e9),
        (Mars, 9.78e8),
        (Jupiter, 5.04e10),
        (Saturn, 6.15e10),
        (Uranus, 6.67e10),
        (Neptune, 1.15e11),
    ],
)
def test_hill_radius(body, expected_r_SOI):
    if expected_r_SOI is not None:
        expected_r_SOI = expected_r_SOI * u.m

    r_SOI = hill_radius(body)

    assert_quantity_allclose(r_SOI, expected_r_SOI, rtol=1e-1)


@pytest.mark.parametrize("missing_body", [Moon, Pluto])
def test_compute_missing_body_hill_radius_raises_error(missing_body):
    with pytest.raises(RuntimeError) as excinfo:
        hill_radius(missing_body)
    assert (
        "To compute the semimajor axis or eccentricity for Moon and Pluto use the JPL ephemeris"
        in excinfo.exconly()
    )


def test_hill_radius_given_a():
    parent = Body(None, 1 * u.km ** 3 / u.s ** 2, "Parent")
    body = Body(parent, 1 * u.km ** 3 / u.s ** 2, "Body")
    r_SOI = hill_radius(body, 1 * u.km, 0.25 * u.one)
    expected_r_SOI = 520.02096 * u.m
    assert_quantity_allclose(r_SOI, expected_r_SOI, rtol=1e-8)
