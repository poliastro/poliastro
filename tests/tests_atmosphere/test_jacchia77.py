import pytest
import numpy as np
from astropy.tests.helper import assert_quantity_allclose

from poliastro.atmosphere.jacchia import Jacchia77

# SOLUTIONS DIRECTLY TAKEN FROM JACCHIA77 IMPLEMENTATION
jacchia77_solutions = {
    0: [288.15, 25.2987, 24.7273, -9.9000, 23.3765, 20.1255, -9.9000, 25.4060, 28.960],
    10: [223.25, 24.8271, 24.2557, -9.9000, 22.9049, 19.6539, -9.9000, 24.9344, 28.960],
    20: [216.65, 24.1595, 23.5881, -9.9000, 22.2373, 18.9863, -9.9000, 24.2668, 28.960],
    30: [226.51, 23.4756, 22.9042, -9.9000, 21.5534, 18.3024, -9.9000, 23.5829, 28.960],
    40: [250.35, 22.8122, 22.2407, -2.7546, 20.8899, 17.6390, -9.9000, 22.9195, 28.960],
    50: [270.65, 22.2221, 21.6507,  3.6757, 20.2999, 17.0489, -9.9000, 22.3294, 28.960],
    60: [247.02, 21.7015, 21.1300,  8.9864, 19.7792, 16.5282, -9.9000, 21.8087, 28.960],
    70: [219.58, 21.1287, 20.5573, 13.0560, 19.2065, 15.9555, -9.9000, 21.2360, 28.960],
    80: [198.63, 20.4767, 19.9052, 15.8573, 18.5544, 15.3034, -9.9000, 20.5839, 28.960],
    85: [188.89, 20.1253, 19.5534, 16.7867, 18.2031, 14.9521, -9.9000, 20.2326, 28.955],
    90: [188.00, 19.7457, 19.1692, 17.3899, 17.8235, 14.5725, -9.9000, 19.8534, 28.908],
    100: [193.69, 18.9709, 18.3232, 17.6649, 17.0487, 13.7977, -9.9000, 19.0803, 28.360],
    120: [350.50, 17.5778, 16.6364, 16.9849, 15.1718, 13.4758, -9.9000, 17.7157, 26.145],
    150: [669.79, 16.4756, 15.4121, 16.2375, 13.7198, 13.1841, 11.7563, 16.6974, 24.059],
    200: [884.39, 15.5002, 14.3152, 15.6287, 12.3803, 12.9872, 11.3932, 15.8828, 21.402],
    300: [975.08, 14.0483, 12.6628, 14.7813, 10.3280, 12.7595, 11.1996, 14.8613, 17.851],
    400: [991.74, 12.7325, 11.1608, 14.0266,  8.4547, 12.5680, 11.1295, 14.0633, 16.179],
    500: [996.42, 11.4719,  9.7212, 13.3058,  6.6580, 12.3869, 11.0792, 13.3632, 14.813],
    560: [997.64, 10.7357,  8.8803, 12.8850,  5.6083, 12.2814, 11.0517, 12.9892, 13.543],
    800: [999.33,  7.9200,  5.6641, 11.2766,  1.5933, 11.8788, 10.9493, 12.0147,  5.939],
    1000: [999.69, 5.7162, 3.1468, 10.0178, -1.5493, 11.5638, 10.8698, 11.6540,  3.788],
    1500: [999.92,  0.6993, -2.5839,  7.1525,  -8.7035, 10.8469, 10.6890, 11.0762,  2.776],
    2000: [999.97, -3.7167, -7.6281,  4.6303, -9.9000, 10.2160, 10.5301, 10.7018,  1.986],
    2400: [999.98, -6.8863, -9.9000,  2.8201, -9.9000,  9.7631, 10.4160, 10.5032,  1.553],
}

@pytest.mark.parametrize("z", jacchia77_solutions.keys())
def test_jacchia77(z):
    #  Z, T, CN2, CO2, CO, CAr, CHe, CH, CM, WM)
    expected_T = jacchia77_solutions[z][0]
    expected_CN2 = jacchia77_solutions[z][1]
    expected_CO2 = jacchia77_solutions[z][2]
    expected_CO = jacchia77_solutions[z][3]
    expected_CAr = jacchia77_solutions[z][4]
    expected_CHe = jacchia77_solutions[z][5]
    expected_CH = jacchia77_solutions[z][6]
    expected_CM = jacchia77_solutions[z][7]
    expected_WM = jacchia77_solutions[z][8]

    properties = Jacchia77().jprop(z, 1000)
    # print(properties)
    for i in range(2,len(properties) - 1):
        if properties[i] > 1.26E-16:
            properties[i] = np.log10(properties[i]) + 6
        else:
            properties[i] = -9.9
    # print(properties)
    Z, T, CN2, CO2, CO, CAr, CHe, CH, CM, WM = properties

    assert_quantity_allclose(T, expected_T, rtol=1e-4)
    assert_quantity_allclose(CN2, expected_CN2, rtol=1e-4)
    assert_quantity_allclose(CO2, expected_CO2, rtol=1e-4)
    assert_quantity_allclose(CO, expected_CO, rtol=1e-4)
    assert_quantity_allclose(CAr, expected_CAr, rtol=1e-4)
    assert_quantity_allclose(CHe, expected_CHe, rtol=1e-4)
    assert_quantity_allclose(CH, expected_CH, rtol=1e-4)
    assert_quantity_allclose(CM, expected_CM, rtol=1e-4)
    assert_quantity_allclose(WM, expected_WM, rtol=1e-3)





