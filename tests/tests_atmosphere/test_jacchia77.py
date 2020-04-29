import numpy as np
import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from poliastro.atmosphere.jacchia import Jacchia77

kmol = u.def_unit("kmol", 1000 * u.mol)

# SOLUTIONS DIRECTLY TAKEN FROM JACCHIA77 REPORT AND
# https://ccmc.gsfc.nasa.gov/pub/modelweb/atmospheric/jacchia/jacchia-77/t1000.out
jacchia77_solutions = {
    90
    * u.km: [
        188.0,
        19.746,
        19.170,
        17.390,
        17.824,
        14.573,
        -9.9000,
        19.854,
        28.91,
        -0.732,
        -5.465,
    ],
    100
    * u.km: [
        193.7,
        18.971,
        18.323,
        17.665,
        17.049,
        13.798,
        -9.9000,
        19.081,
        28.36,
        -1.492,
        -6.247,
    ],
    120
    * u.km: [
        350.5,
        17.578,
        16.637,
        16.985,
        15.172,
        13.476,
        -9.9000,
        17.716,
        26.15,
        -2.599,
        -7.646,
    ],
    150
    * u.km: [
        669.8,
        16.476,
        15.412,
        16.238,
        13.720,
        13.184,
        11.756,
        16.698,
        24.06,
        -3.336,
        -8.701,
    ],
    200
    * u.km: [
        884.4,
        15.501,
        14.315,
        15.629,
        12.381,
        12.987,
        11.392,
        15.883,
        21.40,
        -4.030,
        -9.566,
    ],
    300
    * u.km: [
        975.1,
        14.049,
        12.663,
        14.781,
        10.328,
        12.760,
        11.199,
        14.862,
        17.85,
        -5.009,
        -10.667,
    ],
    400
    * u.km: [
        991.7,
        12.733,
        11.161,
        14.027,
        8.455,
        12.568,
        11.129,
        14.064,
        16.18,
        -5.800,
        -11.507,
    ],
    500
    * u.km: [
        996.4,
        11.472,
        9.722,
        13.306,
        6.658,
        12.387,
        11.079,
        13.363,
        14.81,
        -6.498,
        -12.246,
    ],
    560
    * u.km: [
        997.64,
        10.7357,
        8.8803,
        12.8850,
        5.6083,
        12.2814,
        11.0517,
        12.9892,
        13.543,
        -6.872,
        -12.659,
    ],
    800
    * u.km: [
        999.33,
        7.9200,
        5.6641,
        11.2766,
        1.5933,
        11.8788,
        10.9493,
        12.0147,
        5.939,
        -7.845,
        -13.991,
    ],
    1000
    * u.km: [
        999.69,
        5.7162,
        3.1468,
        10.0178,
        -1.5493,
        11.5638,
        10.8698,
        11.6540,
        3.788,
        -8.206,
        -14.547,
    ],
    1500
    * u.km: [
        999.92,
        0.6993,
        -2.5839,
        7.1525,
        -8.7035,
        10.8469,
        10.6890,
        11.0762,
        2.776,
        -8.783,
        -15.260,
    ],
    2000
    * u.km: [
        999.97,
        -3.7167,
        -7.6281,
        4.6303,
        -9.9000,
        10.2160,
        10.5301,
        10.7018,
        1.986,
        -9.158,
        -15.780,
    ],
    2400
    * u.km: [
        999.98,
        -6.8863,
        -9.9000,
        2.8201,
        -9.9000,
        9.7631,
        10.4160,
        10.5032,
        1.553,
        -9.357,
        -16.148,
    ],
}


@pytest.mark.parametrize("z", jacchia77_solutions.keys())
def test_jacchia77(z):
    #  Z, T, CN2, CO2, CO, CAr, CHe, CH, CM, WM)
    expected_CN2 = jacchia77_solutions[z][1] * (u.m) ** -3
    expected_CO2 = jacchia77_solutions[z][2] * (u.m) ** -3
    expected_CO = jacchia77_solutions[z][3] * (u.m) ** -3
    expected_CAr = jacchia77_solutions[z][4] * (u.m) ** -3
    expected_CHe = jacchia77_solutions[z][5] * (u.m) ** -3
    expected_CH = jacchia77_solutions[z][6] * (u.m) ** -3
    expected_CM = jacchia77_solutions[z][7] * (u.m) ** -3
    expected_WM = jacchia77_solutions[z][8] * (u.kg / kmol)

    properties = Jacchia77().altitude_profile(z, 1000 * u.K)

    for i in range(2, len(properties) - 1):
        if properties[i].value > 1.26e-10:
            properties[i] = np.log10(properties[i].value) * properties[i].unit
        else:
            properties[i] = -9.9 * properties[i].unit

    Z, T, CN2, CO2, CO, CAr, CHe, CH, CM, WM = properties

    assert_quantity_allclose(CN2, expected_CN2, rtol=1e-4)
    assert_quantity_allclose(CO2, expected_CO2, rtol=1e-4)
    assert_quantity_allclose(CO, expected_CO, rtol=1e-4)
    assert_quantity_allclose(CAr, expected_CAr, rtol=1e-4)
    assert_quantity_allclose(CHe, expected_CHe, rtol=1e-4)
    assert_quantity_allclose(CH, expected_CH, rtol=1e-3)
    assert_quantity_allclose(CM, expected_CM, rtol=1e-4)
    assert_quantity_allclose(WM, expected_WM, rtol=1e-3)


@pytest.mark.parametrize("z", jacchia77_solutions.keys())
def test_tempertaure(z):
    expected_T = jacchia77_solutions[z][0] * u.K
    T = Jacchia77().temperature(z, 1000 * u.K)

    assert_quantity_allclose(T, expected_T, rtol=1e-4)


@pytest.mark.parametrize("z", jacchia77_solutions.keys())
def test_pressure(z):
    expected_p = jacchia77_solutions[z][9] * (u.N / u.m / u.m)
    pressure = Jacchia77().pressure(z, 1000 * u.K)
    p = np.log10(pressure.value) * pressure.unit

    assert_quantity_allclose(p, expected_p, rtol=1e-3)


@pytest.mark.parametrize("z", jacchia77_solutions.keys())
def test_density(z):
    expected_rho = jacchia77_solutions[z][10] * (u.kg / u.m / u.m / u.m)
    density = Jacchia77().density(z, 1000 * u.K)
    rho = np.log10(density.value) * density.unit

    assert_quantity_allclose(rho, expected_rho, rtol=1e-2)


def test_outside_upper_limit_coesa76():
    with pytest.raises(ValueError) as excinfo:
        alt = 2501.0 * u.km
        Texo = 1000 * u.K
        Jacchia77().altitude_profile(alt, Texo)
    assert (
        "ValueError: Jacchia77 has been implemented in range 90km - 2500km."
        in excinfo.exconly()
    )


def test_outside_lower_limit_coesa76():
    with pytest.raises(ValueError) as excinfo:
        alt = 89.0 * u.km
        Texo = 1000 * u.K
        Jacchia77().altitude_profile(alt, Texo)
    assert (
        "ValueError: Jacchia77 has been implemented in range 90km - 2500km."
        in excinfo.exconly()
    )
