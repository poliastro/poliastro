import numpy as np
from numpy import radians
from numpy.testing import TestCase, assert_array_almost_equal, \
    assert_almost_equal, run_module_suite

from poliastro.constants import k_Earth, k_Sun

from poliastro.iod import lambert


class TestLambert(TestCase):
    def test_vallado75(self):
        k = k_Earth
        r0 = np.array([15945.34, 0.0, 0.0])
        r = np.array([12214.83399, 10249.46731, 0.0])
        tof = 76.0 * 60.0
        va, vb = lambert(k, r0, r, tof)
        assert_array_almost_equal(va, np.array([2.058925, 2.915956, 0.0]), decimal=4)
        assert_array_almost_equal(vb, np.array([-3.451569, 0.910301, 0.0]), decimal=5)

    def test_curtis52(self):
        k = k_Earth
        r0 = np.array([5000.0, 10000.0, 2100.0])
        r = np.array([-14600.0, 2500.0, 7000.0])
        tof = 3600.0
        va, vb = lambert(k, r0, r, tof)
        assert_array_almost_equal(va, np.array([-5.9925, 1.9254, 3.2456]), decimal=4)
        assert_array_almost_equal(vb, np.array([-3.3125, -4.1966, -0.38529]), decimal=4)

    def test_curtis53(self):
        k = k_Earth
        r0 = np.array([273378.0, 0.0, 0.0])
        r = np.array([145820.0, 12758.0, 0.0])
        tof = 13.5 * 3600.0
        va, vb = lambert(k, r0, r, tof)
        # ERRATA: j component is positive
        assert_array_almost_equal(va, np.array([-2.4356, 0.26741, 0.0]), decimal=3)


if __name__ == '__main__':
    run_module_suite()
