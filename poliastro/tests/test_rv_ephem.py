from datetime import datetime

import numpy as np
from numpy import radians
from numpy.testing import TestCase, assert_array_almost_equal, \
    run_module_suite

from poliastro import ephem, twobody
from poliastro.constants import AU, k_Sun


class TestRVEphem(TestCase):
    def test_vallado55(self):
        # Data from HORIZONS, the xy-plane is the ecliptic
        TU = 86400
        jday = ephem.jd(datetime(1994, 5, 20, 20))
        coe = ephem.mean_elements(jday, ephem.JUPITER)
        r_XYZ, v_XYZ = twobody.coe2rv(k_Sun, *coe)
        r0_XYZ = np.array([-4.0686995, -3.5897288, 0.1059762]) * AU
        v0_XYZ = np.array([0.0048974, -0.0053138, -0.00008791]) * AU / TU
        assert_array_almost_equal(r_XYZ / AU, r0_XYZ / AU, decimal=2)
        assert_array_almost_equal(v_XYZ, v0_XYZ, decimal=2)

if __name__ == '__main__':
    run_module_suite()
