import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from poliastro.earth.atmosphere import COESA76
from poliastro.earth.atmosphere.coesa76 import p_coeff, rho_coeff

coesa76 = COESA76()


def test_outside_altitude_range_coesa76():
    with pytest.raises(ValueError) as excinfo:
        r0 = 6356.766 * u.km
        coesa76._check_altitude(1001 * u.km, r0)
    assert (
        "ValueError: Geometric altitude must be in range [0.0 km, 1000.0 km]"
        in excinfo.exconly()
    )


def test_get_index_coesa76():
    expected_i = 7
    z = 86 * u.km
    i = coesa76._get_index(z, coesa76.zb_levels)
    assert i == expected_i


def test_coefficients_over_86km():
    # Expected pressure coefficients
    expected_p = [
        9.814674e-11,
        -1.654439e-07,
        1.148115e-04,
        -0.05431334,
        -2.011365,
    ]
    expected_rho = [
        1.140564e-10,
        -2.130756e-07,
        1.570762e-04,
        -0.07029296,
        -12.89844,
    ]

    assert (
        coesa76._get_coefficients_avobe_86(350 * u.km, p_coeff) == expected_p
    )
    assert (
        coesa76._get_coefficients_avobe_86(350 * u.km, rho_coeff)
        == expected_rho
    )


# SOLUTIONS DIRECTLY TAKEN FROM COESA76 REPORT
coesa76_solutions = {
    0.5 * u.km: [284.90 * u.K, 9.5461e2 * u.mbar, 1.1673 * u.kg / u.m**3],
    1.0 * u.km: [281.651 * u.K, 8.9876e2 * u.mbar, 1.1117 * u.kg / u.m**3],
    10 * u.km: [223.252 * u.K, 2.6499e2 * u.mbar, 4.1351e-1 * u.kg / u.m**3],
    77
    * u.km: [204.493 * u.K, 1.7286e-2 * u.mbar, 2.9448e-5 * u.kg / u.m**3],
    86 * u.km: [186.87 * u.K, 3.7338e-3 * u.mbar, 6.958e-6 * u.kg / u.m**3],
    92 * u.km: [186.96 * u.K, 1.2887e-3 * u.mbar, 2.393e-6 * u.kg / u.m**3],
    230
    * u.km: [915.78 * u.K, 3.9276e-7 * u.mbar, 1.029e-10 * u.kg / u.m**3],
    1000
    * u.km: [1000.0 * u.K, 7.5138e-11 * u.mbar, 3.561e-15 * u.kg / u.m**3],
}


@pytest.mark.parametrize("z", coesa76_solutions.keys())
def test_properties_coesa76(z):
    # Get expected values from official data
    expected_T = coesa76_solutions[z][0]
    expected_p = coesa76_solutions[z][1]
    expected_rho = coesa76_solutions[z][2]

    T, p, rho = coesa76.properties(z)

    assert_quantity_allclose(T, expected_T, rtol=1e-4)
    assert_quantity_allclose(p, expected_p, rtol=1e-4)
    assert_quantity_allclose(rho, expected_rho, rtol=1e-3)


# DATA DIRECTLY TAKEN FROM TABLE-III COESA76 REPORT
sound_speed_viscosity_conductivity = {
    0.5
    * u.km: [
        338.37 * (u.m / u.s),
        1.7737e-5 * (u.N * u.s / (u.m) ** 2),
        2.5106e-2 * (u.J / u.m / u.s / u.K),
    ],
    10
    * u.km: [
        299.53 * (u.m / u.s),
        1.4577e-5 * (u.N * u.s / (u.m) ** 2),
        2.0088e-2 * (u.J / u.m / u.s / u.K),
    ],
    24
    * u.km: [
        297.72 * (u.m / u.s),
        1.4430e-5 * (u.N * u.s / (u.m) ** 2),
        1.9862e-2 * (u.J / u.m / u.s / u.K),
    ],
    41
    * u.km: [
        318.94 * (u.m / u.s),
        1.6151e-5 * (u.N * u.s / (u.m) ** 2),
        2.2556e-2 * (u.J / u.m / u.s / u.K),
    ],
    50
    * u.km: [
        329.80 * (u.m / u.s),
        1.7037e-5 * (u.N * u.s / (u.m) ** 2),
        2.3973e-2 * (u.J / u.m / u.s / u.K),
    ],
    67
    * u.km: [
        302.57 * (u.m / u.s),
        1.4823e-5 * (u.N * u.s / (u.m) ** 2),
        2.0469e-2 * (u.J / u.m / u.s / u.K),
    ],
    85
    * u.km: [
        275.52 * (u.m / u.s),
        1.2647e-5 * (u.N * u.s / (u.m) ** 2),
        1.7162e-2 * (u.J / u.m / u.s / u.K),
    ],
}


@pytest.mark.parametrize("z", sound_speed_viscosity_conductivity.keys())
def test_sound_speed_viscosity_conductivity(z):
    expected_Cs = sound_speed_viscosity_conductivity[z][0]
    expected_mu = sound_speed_viscosity_conductivity[z][1]
    expected_k = sound_speed_viscosity_conductivity[z][2]

    Cs = coesa76.sound_speed(z)
    mu = coesa76.viscosity(z)
    k = coesa76.thermal_conductivity(z)

    assert_quantity_allclose(Cs, expected_Cs, rtol=1e-4)
    assert_quantity_allclose(mu, expected_mu, rtol=1e-4)
    assert_quantity_allclose(k, expected_k, rtol=1e-2)


def test_sound_speed_over_86km():
    z = 87 * u.km
    # Test speed of sound over 86 km
    with pytest.raises(ValueError) as excinfo:
        coesa76.sound_speed(z)
    assert (
        "ValueError: Speed of sound in COESA76 has just been implemented up to 86km."
        in excinfo.exconly()
    )


def test_viscosity_over_86km():
    z = 87 * u.km
    # Test viscosity over 86 km
    with pytest.raises(ValueError) as excinfo:
        coesa76.viscosity(z)
    assert (
        "ValueError: Dynamic Viscosity in COESA76 has just been implemented up to 86km."
        in excinfo.exconly()
    )


def test_conductivity_over_86km():
    z = 87 * u.km
    # Test thermal conductivity over 86 km
    with pytest.raises(ValueError) as excinfo:
        coesa76.thermal_conductivity(z)
    assert (
        "ValueError: Thermal conductivity in COESA76 has just been implemented up to 86km."
        in excinfo.exconly()
    )
