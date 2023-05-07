from math import cos, pi, sin

from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from poliastro.bodies import Earth, Moon
from poliastro.threebody.restricted import lagrange_points_vec


def test_lagrange_points_vec():
    # Figure 2.36 from Curtis

    deg60 = 60 * pi / 180
    expected_L1 = 326400 * ([1, 0, 0] * u.km)
    expected_L2 = 449100 * ([1, 0, 0] * u.km)
    expected_L3 = -381600 * ([1, 0, 0] * u.km)
    expected_L4 = 384400 * ([cos(deg60), sin(deg60), 0] * u.km)
    expected_L5 = 384400 * ([cos(deg60), -sin(deg60), 0] * u.km)

    earth_mass = Earth.mass
    moon_mass = Moon.mass

    # Values from Curtis
    # earth_mass = 5.974e24 * u.kg
    # moon_mass = 73.48e21 * u.kg

    L1, L2, L3, L4, L5 = lagrange_points_vec(
        m1=earth_mass,
        r1=([0, 0, 0] * u.km),
        m2=moon_mass,
        r2=384400 * ([1, 0, 0] * u.km),
        n=[0, 0, 1] * u.one,
    )

    assert_quantity_allclose(L1, expected_L1, rtol=1.0e-3)
    assert_quantity_allclose(L2, expected_L2, rtol=1.0e-3)
    assert_quantity_allclose(L3, expected_L3, rtol=1.0e-3)
    assert_quantity_allclose(L4, expected_L4, rtol=1.0e-3)
    assert_quantity_allclose(L5, expected_L5, rtol=1.0e-3)
