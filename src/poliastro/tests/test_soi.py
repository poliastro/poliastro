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
from poliastro.threebody.soi import laplace_radius


def test_laplace_radius():
    # Data from Table A.2., Curtis "Orbital Mechanics for Engineering Students"
    data = [
        # body, SOI radius (m)
        (Sun, None),
        (Mercury, 1.12e8),
        (Venus, 6.16e8),
        (Earth, 9.25e8),
        # (Moon, 6.61e7),
        (Mars, 5.77e8),
        (Jupiter, 4.82e10),
        (Saturn, 5.48e10),
        (Uranus, 5.18e10),
        (Neptune, 8.66e10),
        # (Pluto, 3.08e9)
    ]
    for row in data:

        body, expected_r_SOI = row
        if expected_r_SOI is not None:
            expected_r_SOI = expected_r_SOI * u.m
        else:
            continue

        r_SOI = laplace_radius(body)

        assert_quantity_allclose(r_SOI, expected_r_SOI, rtol=1e-1)


"""def test_hill_radius():
    # Data from Table A.2., Curtis "Orbital Mechanics for Engineering Students"
    data = [
        # body, SOI radius (m)
        (Sun, None),
        (Mercury, 1.12e8),
        (Venus, 6.16e8),
        (Earth, 9.25e8),
        # (Moon, 6.61e7),
        (Mars, 5.77e8),
        (Jupiter, 4.82e10),
        (Saturn, 5.48e10),
        (Uranus, 5.18e10),
        (Neptune, 8.66e10),
        # (Pluto, 3.08e9)
    ]
    for row in data:

        body, expected_r_HR = row
        if expected_r_HR is not None:
            expected_r_HR = expected_r_HR * u.m
        else:
            continue

        r_HR = hill_radius(body)

        assert_quantity_allclose(r_SOI, expected_r_SOI, rtol=1e-1)
"""


@pytest.mark.parametrize("missing_body", [Moon, Pluto])
def test_compute_missing_body_soi_raises_error(missing_body):
    with pytest.raises(RuntimeError) as excinfo:
        laplace_radius(missing_body)
    assert (
        "To compute the semimajor axis for Moon and Pluto use the JPL ephemeris"
        in excinfo.exconly()
    )


def test_laplace_radius_given_a():
    parent = Body(None, 1 * u.km ** 3 / u.s ** 2, "Parent")
    body = Body(parent, 1 * u.km ** 3 / u.s ** 2, "Body")
    r_SOI = laplace_radius(body, 1 * u.km)

    assert r_SOI == 1 * u.km


def test_hill_radius_given_a():
    parent = Body(None, 1 * u.km ** 3 / u.s ** 2, "Parent")
    body = Body(parent, 1 * u.km ** 3 / u.s ** 2, "Body")
    r_HR = hill_radius(body, 1 * u.km)

    assert r_HR == ((1 / 3) * u.km) ** (1 / 3)
