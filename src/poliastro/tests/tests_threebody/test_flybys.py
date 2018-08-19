import pytest

from astropy.tests.helper import assert_quantity_allclose
from astropy import units as u

from poliastro.bodies import Venus
from poliastro.threebody.flybys import compute_flyby


@pytest.mark.parametrize("theta, expected_V_2_v", [
    (0 * u.deg, [31.73, 1.766, 0] * u.km / u.s),
    (180 * u.deg, [37.14, -3.074, 0] * u.km / u.s),
])
def test_flyby_curtis(theta, expected_V_2_v):
    # Data from Curtis, Example 8.6
    periapsis_h = 300 * u.km
    V_1_v = [37.51, 2.782, 0] * u.km / u.s  # In some reference frame
    V = [35.02, 0, 0] * u.km / u.s

    expected_delta = 103.6 * u.deg

    V_2_v, delta = compute_flyby(V_1_v, V, Venus.k, Venus.R + periapsis_h, theta)

    assert_quantity_allclose(V_2_v, expected_V_2_v, rtol=1e-3, atol=1e-15 * u.km / u.s)
    assert_quantity_allclose(delta, expected_delta, rtol=1e-3)
