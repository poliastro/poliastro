from datetime import datetime

import numpy as np
from numpy import radians
from numpy.testing import TestCase, assert_almost_equal, run_module_suite

from poliastro import ephem
from poliastro.constants import AU


class TestEphem(TestCase):
    def test_vallado55(self):
        dd = datetime(1994, 5, 20, 20)
        a, ecc, inc, omega, argp, nu = ephem.mean_elements(dd, ephem.JUPITER)
        assert_almost_equal(a / AU, 5.202603, decimal=6)
        assert_almost_equal(ecc, 0.048486, decimal=6)
        assert_almost_equal(inc, radians(1.303270), decimal=4)
        assert_almost_equal(omega, radians(100.454519), decimal=2)
        assert_almost_equal(argp, radians(-86.135316), decimal=2)
        assert_almost_equal(nu, radians(206.95453))


class TestJd(TestCase):
    def test_vallado55(self):
        dd = datetime(1994, 5, 20, 20)
        assert_almost_equal(ephem.jd(dd), 2449493.333, decimal=3)


if __name__ == '__main__':
    run_module_suite()
