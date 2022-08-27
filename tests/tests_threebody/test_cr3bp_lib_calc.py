import numpy as np
import pytest

# from astropy import units as u
# from astropy.units import L_ND
from astropy.tests.helper import assert_quantity_allclose

from poliastro.threebody.cr3bp_char_quant import SystemChars
from poliastro.threebody.cr3bp_lib_calc import lib_pt_loc


@pytest.mark.parametrize(
    "SysChars, conv_tol, expected_lib_pt_loc",
    [
        (
            SystemChars(
                "EarthMoon",
                1.215058560962404e-02,
                389703.2648292776,
                382981.2891290545,
            ),
            1e-12,
            np.array(
                [
                    [0.8369151257723572, 0.0, 0.0],
                    [1.155682165444884, 0.0, 0.0],
                    [-1.005062645810278, 0.0, 0.0],
                    [0.4878494143903759, 0.8660254037844386, 0.0],
                    [0.4878494143903759, -0.8660254037844386, 0.0],
                ]
            ),
        ),
        # Earth-Moon mu, l*, t*, libration ponints: https://ssd-api.jpl.nasa.gov/periodic_orbits.api?sys=earth-moon&family=halo&libr=1&branch=S
    ],
)
def test_lib_pt_loc(SysChars, conv_tol, expected_lib_pt_loc):

    lib_pt = lib_pt_loc(SysChars, conv_tol)

    # L_ND = u.def_unit("dist_nd", SysChars.lstar)

    lib_pt = lib_pt.value

    assert_quantity_allclose(lib_pt, expected_lib_pt_loc, 1e-6)
