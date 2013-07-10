import numpy as np
from numpy import radians
from numpy.testing import TestCase, assert_array_almost_equal, \
    assert_almost_equal, run_module_suite

from poliastro.angles import M2nu, nu2M


class TestM2nu(TestCase):
    def test_vectorize(self):
        N = 10
        M_vec = np.linspace(0, 2 * np.pi, N)
        ecc = 0.1
        e0, nu = M2nu(ecc, M_vec)
        assert e0.shape == M_vec.shape
        assert nu.shape == M_vec.shape

    def test_data(self):
        # Data from Schlesinger & Udick, 1912
        data = [
            # ecc, M (deg), nu (deg)
            (0.0, 0.0, 0.0),
            (0.05, 10.0, 11.06),
            (0.06, 30.0, 33.67),
            (0.04, 120.0, 123.87),
            (0.14, 65.0, 80.50),
            (0.19, 21.0, 30.94),
            (0.35, 65.0, 105.71),
            (0.48, 180.0, 180.0),
            (0.75, 125.0, 167.57)
        ]
        for row in data:
            ecc, M0, nu0 = row
            _, nu = M2nu(ecc, radians(M0))
            assert_almost_equal(nu, radians(nu0), decimal=3)


class TestNu2m(TestCase):
    def test_vectorize(self):
        N = 10
        nu_vec = np.linspace(0, 2 * np.pi, N)
        ecc = 0.1
        e0, M = nu2M(ecc, nu_vec)
        assert e0.shape == nu_vec.shape
        assert M.shape == nu_vec.shape


if __name__ == '__main__':
    run_module_suite()
