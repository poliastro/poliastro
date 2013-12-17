import numpy as np
from numpy import radians
from numpy.testing import TestCase, assert_array_almost_equal, \
    assert_almost_equal, run_module_suite

from astropy import units
from astropy.constants import G, M_earth, R_earth

from poliastro.twobody import hohmann, bielliptic, coe2rv, rv2coe, kepler

# TODO: Low precision in some results, why? Use canonical units?
# TODO: Test for exceptions

k_earth = (G * M_earth)


class TestHohmann(TestCase):
    def test_vallado61(self):
        alt_i = 191.34411  # km
        alt_f = 35781.34857  # km
        k = k_earth.to(units.km ** 3 / units.s ** 2).value
        R = R_earth.to(units.km).value
        expected_dv = 3.935224  # km/s
        expected_t_trans = 5.256713  # h
        r_i = R + alt_i
        r_f = R + alt_f
        dva, dvb, _, t_trans = hohmann(k, r_i, r_f)
        dv = abs(dva) + abs(dvb)
        assert_almost_equal(dv, expected_dv, decimal=6)
        assert_almost_equal(t_trans / 3600, expected_t_trans, decimal=6)


class TestBielliptic(TestCase):
    def test_vallado62(self):
        alt_i = 191.34411  # km
        alt_b = 503873  # km
        alt_f = 376310.0  # km
        k = k_earth.to(units.km ** 3 / units.s ** 2).value
        R = R_earth.to(units.km).value
        expected_dv = 3.904057  # km/s
        expected_t_trans = 593.919803  # h
        r_i = R + alt_i
        r_b = R + alt_b
        r_f = R + alt_f
        dva, dvb, dvc, _, _, t_trans1, t_trans2 = bielliptic(k, r_i, r_b, r_f)
        dv = abs(dva) + abs(dvb) + abs(dvc)
        t_trans = t_trans1 + t_trans2
        assert_almost_equal(dv, expected_dv, decimal=6)
        assert_almost_equal(t_trans / 3600, expected_t_trans, decimal=3)


class TestCoe2rv(TestCase):
    def test_vectorize(self):
        N = 50
        k = k_earth.to(units.km ** 3 / units.s ** 2).value
        p = 11067.790
        ecc = 0.83285
        a = p / (1 - ecc ** 2)
        inc = radians(87.87)
        omega = radians(227.89)
        argp = radians(53.38)
        nu = np.linspace(0, 2 * np.pi, num=N)
        r, v = coe2rv(k, a, ecc, inc, omega, argp, nu)
        assert r.shape, (3,) == nu.shape
        assert v.shape, (3,) == nu.shape

    def test_vallado26(self):
        k = k_earth.to(units.km ** 3 / units.s ** 2).value
        p = 11067.790
        ecc = 0.83285
        a = p / (1 - ecc ** 2)
        inc = radians(87.87)
        omega = radians(227.89)
        argp = radians(53.38)
        nu = radians(92.335)
        r, v = coe2rv(k, a, ecc, inc, omega, argp, nu)
        assert_array_almost_equal(r, np.array([6525.344, 6861.535, 6449.125]),
                                  decimal=1)
        assert_array_almost_equal(v, np.array([4.902276, 5.533124, -1.975709]),
                                  decimal=4)


class TestRv2coe(TestCase):
    def test_vallado25(self):
        k = k_earth.to(units.km ** 3 / units.s ** 2).value
        r = np.array([6524.384, 6862.875, 6448.296])
        v = np.array([4.901327, 5.533756, -1.976341])
        a, ecc, inc, omega, argp, nu = rv2coe(k, r, v)
        p = a * (1 - ecc ** 2)
        assert_almost_equal(p, 11067.79, decimal=0)
        assert_almost_equal(ecc, 0.832853, decimal=4)
        assert_almost_equal(inc, radians(87.870), decimal=4)
        assert_almost_equal(omega, radians(227.89), decimal=3)
        assert_almost_equal(argp, radians(53.38), decimal=3)
        assert_almost_equal(nu, radians(92.335), decimal=5)


class TestKepler(TestCase):
    def test_vallado24(self):
        k = k_earth.to(units.km ** 3 / units.s ** 2).value
        r0 = np.array([1131.340, -2282.343, 6672.423])
        v0 = np.array([-5.64305, 4.30333, 2.42879])
        tof = 40 * 60.0
        r, v = kepler(k, r0, v0, tof)
        assert_array_almost_equal(r, np.array([-4219.7527, 4363.0292, -3958.7666]),
                                  decimal=4)
        assert_array_almost_equal(v, np.array([3.689866, -1.916735, -6.112511]))


if __name__ == '__main__':
    run_module_suite()
