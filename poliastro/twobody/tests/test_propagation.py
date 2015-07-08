# coding: utf-8
import numpy as np
from numpy.testing import assert_array_almost_equal

from astropy import units as u

from poliastro.bodies import Earth

from poliastro.twobody.propagation import cowell


def test_cowell_propagation_with_zero_acceleration_equals_kepler():
    # Data from Vallado, example 2.4
    k = Earth.k.to(u.km**3 / u.s**2).value

    r0 = np.array([1131.340, -2282.343, 6672.423])  # km
    v0 = np.array([-5.64305, 4.30333, 2.42879])  # km/s
    tof = 40 * 60.0  # s

    expected_r = np.array([-4219.7527, 4363.0292, -3958.7666])
    expected_v = np.array([3.689866, -1.916735, -6.112511])

    r, v = cowell(k, r0, v0, tof, None)

    assert_array_almost_equal(r, expected_r, decimal=1)
    assert_array_almost_equal(v, expected_v, decimal=4)
