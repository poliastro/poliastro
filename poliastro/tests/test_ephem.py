from datetime import datetime

import numpy as np
from numpy import radians
from numpy.testing import TestCase, assert_almost_equal, run_module_suite

from astropy import units
from astropy.constants import au

from poliastro import ephem

AU = au.to(units.km).value


class TestEphem(TestCase):
    def test_vallado55(self):
        dd = datetime(1994, 5, 20, 20)
        a, ecc, inc, omega, argp, nu = ephem.mean_elements(dd, ephem.JUPITER)
        assert_almost_equal(a / AU, 5.202895, decimal=3)
        assert_almost_equal(ecc, 0.048319, decimal=3)
        assert_almost_equal(inc, 0.022770, decimal=4)
        assert_almost_equal(omega, 1.753543, decimal=2)
        assert_almost_equal(argp, 4.803144, decimal=1)
        assert_almost_equal(nu, 3.590915, decimal=1)


if __name__ == '__main__':
    run_module_suite()
