import numpy as np
from numpy import radians
from numpy.testing import TestCase, assert_array_almost_equal, \
    assert_almost_equal, run_module_suite

from poliastro.constants import k_Earth, k_Sun

from poliastro.twobody import coe2rv, rv2coe, kepler

# TODO: Low precision in some results, why? Use canonical units?


class TestCoe2rv(TestCase):
    def test_vallado26(self):
        k = k_Earth
        p = 11067.790
        ecc = 0.83285
        inc = radians(87.87)
        omega = radians(227.89)
        argp = radians(53.38)
        nu = radians(92.335)
        r, v = coe2rv(k, p, ecc, inc, omega, argp, nu)
        assert_array_almost_equal(r, np.array([6525.344, 6861.535, 6449.125]), decimal=1)
        assert_array_almost_equal(v, np.array([4.902276, 5.533124, -1.975709]), decimal=4)


class TestRv2coe(TestCase):
    def test_vallado25(self):
        k = k_Earth
        r = np.array([6524.384, 6862.875, 6448.296])
        v = np.array([4.901327, 5.533756, -1.976341])
        p, ecc, inc, omega, argp, nu = rv2coe(k, r, v)
        assert_almost_equal(p, 11067.79, decimal=0)
        assert_almost_equal(ecc, 0.832853, decimal=4)
        assert_almost_equal(inc, radians(87.870), decimal=4)
        assert_almost_equal(omega, radians(227.89), decimal=3)
        assert_almost_equal(argp, radians(53.38), decimal=3)
        assert_almost_equal(nu, radians(92.335), decimal=5)


class TestKepler(TestCase):
    def test_vallado24(self):
        k = k_Earth
        r0 = np.array([1131.340, -2282.343, 6672.423])
        v0 = np.array([-5.64305, 4.30333, 2.42879])
        tof = 40 * 60.0
        r, v, error = kepler(k, r0, v0, tof)
        assert_array_almost_equal(r, np.array([-4219.7527, 4363.0292, -3958.7666]), decimal=4)
        assert_array_almost_equal(v, np.array([3.689866, -1.916735, -6.112511]))


if __name__ == '__main__':
    run_module_suite()
