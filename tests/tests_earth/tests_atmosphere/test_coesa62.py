import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from astropy.units import imperial

from poliastro.earth.atmosphere import COESA62

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
    0.5 * u.km: [284.900 * u.K, 9.54612e2 * u.mbar, 1.1673 * u.kg / u.m**3],
    1.0 * u.km: [281.651 * u.K, 8.98762e2 * u.mbar, 1.1117 * u.kg / u.m**3],
    10.0
    * u.km: [223.252 * u.K, 2.64999e2 * u.mbar, 4.1351e-1 * u.kg / u.m**3],
    77.0
    * u.km: [192.340 * u.K, 1.7725e-2 * u.mbar, 3.210e-5 * u.kg / u.m**3],
    86.0
    * u.km: [180.65 * u.K, 3.4313e-3 * u.mbar, 6.617e-6 * u.kg / u.m**3],
    97.0
    * u.km: [201.65 * u.K, 4.8709e-4 * u.mbar, 8.415e-7 * u.kg / u.m**3],
    103.0
    * u.km: [225.65 * u.K, 1.9074e-4 * u.mbar, 2.945e-7 * u.kg / u.m**3],
    115.0
    * u.km: [310.65 * u.K, 4.1224e-5 * u.mbar, 4.623e-8 * u.kg / u.m**3],
    132.0
    * u.km: [600.65 * u.K, 1.0909e-5 * u.mbar, 6.327e-9 * u.kg / u.m**3],
    157.0
    * u.km: [1065.65 * u.K, 4.0409e-6 * u.mbar, 1.321e-9 * u.kg / u.m**3],
    183.0
    * u.km: [1301.65 * u.K, 1.9979e-6 * u.mbar, 5.347e-10 * u.kg / u.m**3],
    201.0
    * u.km: [1405.65 * u.K, 1.3037e-6 * u.mbar, 3.231e-10 * u.kg / u.m**3],
    258.0
    * u.km: [1662.65 * u.K, 4.0061e-7 * u.mbar, 8.394e-11 * u.kg / u.m**3],
    340.0
    * u.km: [1962.65 * u.K, 9.8014e-8 * u.mbar, 1.740e-11 * u.kg / u.m**3],
    482.0
    * u.km: [2373.85 * u.K, 1.3667e-8 * u.mbar, 2.006e-12 * u.kg / u.m**3],
    576.0
    * u.km: [2549.85 * u.K, 4.5072e-9 * u.mbar, 6.158e-13 * u.kg / u.m**3],
    698.0
    * u.km: [2698.45 * u.K, 1.2165e-9 * u.mbar, 1.570e-13 * u.kg / u.m**3],
}


@pytest.mark.parametrize("z", coesa62_solutions.keys())
def test_properties_coesa62(z):
    # Get expected values from official data
    expected_T = coesa62_solutions[z][0]
    expected_p = coesa62_solutions[z][1]
    expected_rho = coesa62_solutions[z][2]

    T, p, rho = coesa62.properties(z)

    assert_quantity_allclose(T, expected_T, rtol=1e-4)
    assert_quantity_allclose(p, expected_p, rtol=1e-3)
    assert_quantity_allclose(rho, expected_rho, rtol=1e-3)


# DATA DIRECTLY TAKEN FROM TABLE-III COESA62 REPORT
sound_speed_viscosity_conductivity = {
    0.5
    * u.km: [
        338.37 * (u.m / u.s),
        1.7737e-5 * (u.kg / u.m / u.s),
        5.9919e-6 * (imperial.kcal / u.m / u.s / u.K),
    ],
    10
    * u.km: [
        299.532 * (u.m / u.s),
        1.4577e-5 * (u.kg / u.m / u.s),
        4.7942e-6 * (imperial.kcal / u.m / u.s / u.K),
    ],
    24
    * u.km: [
        297.72 * (u.m / u.s),
        1.4430e-5 * (u.kg / u.m / u.s),
        4.7403e-6 * (imperial.kcal / u.m / u.s / u.K),
    ],
    41
    * u.km: [
        318.936 * (u.m / u.s),
        1.6151e-5 * (u.kg / u.m / u.s),
        5.3833e-6 * (imperial.kcal / u.m / u.s / u.K),
    ],
    50
    * u.km: [
        329.799 * (u.m / u.s),
        1.7037e-5 * (u.kg / u.m / u.s),
        5.7214e-6 * (imperial.kcal / u.m / u.s / u.K),
    ],
    67
    * u.km: [
        304.979 * (u.m / u.s),
        1.5018e-5 * (u.kg / u.m / u.s),
        4.9575e-6 * (imperial.kcal / u.m / u.s / u.K),
    ],
    89
    * u.km: [
        269.44 * (u.m / u.s),
        1.216e-5 * (u.kg / u.m / u.s),
        3.925e-6 * (imperial.kcal / u.m / u.s / u.K),
    ],
}


@pytest.mark.parametrize("z", sound_speed_viscosity_conductivity.keys())
def test_sound_speed_viscosity_conductivity(z):
    expected_Cs = sound_speed_viscosity_conductivity[z][0]
    expected_mu = sound_speed_viscosity_conductivity[z][1]
    expected_k = sound_speed_viscosity_conductivity[z][2]

    Cs = coesa62.sound_speed(z)
    mu = coesa62.viscosity(z)
    k = coesa62.thermal_conductivity(z)

    assert_quantity_allclose(Cs, expected_Cs, rtol=1e-4)
    assert_quantity_allclose(mu, expected_mu, rtol=1e-3)
    assert_quantity_allclose(k, expected_k, rtol=1e-4)


def test_sound_speed_over_90km():
    z = 91 * u.km
    # Test speed of sound over 90 km
    with pytest.raises(ValueError) as excinfo:
        coesa62.sound_speed(z)
    assert (
        "ValueError: Speed of sound in COESA62 has just been implemented up to 90km."
        in excinfo.exconly()
    )


def test_viscosity_over_90km():
    z = 91 * u.km
    # Test viscosity over 90 km
    with pytest.raises(ValueError) as excinfo:
        coesa62.viscosity(z)
    assert (
        "ValueError: Dynamic Viscosity in COESA62 has just been implemented up to 90km."
        in excinfo.exconly()
    )


def test_conductivity_over_90km():
    z = 91 * u.km
    # Test thermal conductivity over 90 km
    with pytest.raises(ValueError) as excinfo:
        coesa62.thermal_conductivity(z)
    assert (
        "ValueError: Thermal conductivity in COESA62 has just been implemented up to 90km."
        in excinfo.exconly()
    )
