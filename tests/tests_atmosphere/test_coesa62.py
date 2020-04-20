import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from poliastro.atmosphere import COESA62

coesa62 = COESA62()


def test_outside_altitude_range_coesa62():
    with pytest.raises(ValueError) as excinfo:
        r0 = 6356.766 * u.km
        coesa62._check_altitude(701 * u.km, r0)
    assert (
        "ValueError: Geometric altitude must be in range [0.0 km, 700.0 km]"
        in excinfo.exconly()
    )


def test_get_index_coesa62():
    expected_i = 8
    z = 90 * u.km
    i = coesa62._get_index(z, coesa62.zb_levels)
    assert i == expected_i


# SOLUTIONS DIRECTLY TAKEN FROM COESA62 REPORT
coesa62_solutions = {
    0.5 * u.km: [284.900 * u.K, 9.54612e2 * u.mbar, 1.1673 * u.kg / u.m ** 3],
    1.0 * u.km: [281.651 * u.K, 8.98762e2 * u.mbar, 1.1117 * u.kg / u.m ** 3],
    10.0 * u.km: [223.252 * u.K, 2.64999e2 * u.mbar, 4.1351e-1 * u.kg / u.m ** 3],
    77.0 * u.km: [192.340 * u.K, 1.7725e-2 * u.mbar, 3.210e-5 * u.kg / u.m ** 3],
    86.0 * u.km: [180.65 * u.K, 3.4313e-3 * u.mbar, 6.617e-6 * u.kg / u.m ** 3],
}


@pytest.mark.parametrize("z", coesa62_solutions.keys())
def test_properties_coesa62(z):
    # Get expected values from official data
    expected_T = coesa62_solutions[z][0]
    expected_p = coesa62_solutions[z][1]
    expected_rho = coesa62_solutions[z][2]

    T, p, rho = coesa62.properties(z)

    assert_quantity_allclose(T, expected_T, rtol=1e-4)
    assert_quantity_allclose(p, expected_p, rtol=1e-4)
    assert_quantity_allclose(rho, expected_rho, rtol=1e-3)
