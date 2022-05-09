from numpy import cos, cosh, sin, sinh
from numpy.testing import assert_allclose, assert_equal

from poliastro.core.stumpff import c2, c3


def test_stumpff_functions_near_zero():
    psi = 0.5
    expected_c2 = (1 - cos(psi**0.5)) / psi
    expected_c3 = (psi**0.5 - sin(psi**0.5)) / psi**1.5

    assert_allclose(c2(psi), expected_c2)
    assert_allclose(c3(psi), expected_c3)


def test_stumpff_functions_above_zero():
    psi = 3.0
    expected_c2 = (1 - cos(psi**0.5)) / psi
    expected_c3 = (psi**0.5 - sin(psi**0.5)) / psi**1.5

    assert_equal(c2(psi), expected_c2)
    assert_equal(c3(psi), expected_c3)


def test_stumpff_functions_under_zero():
    psi = -3.0
    expected_c2 = (cosh((-psi) ** 0.5) - 1) / (-psi)
    expected_c3 = (sinh((-psi) ** 0.5) - (-psi) ** 0.5) / (-psi) ** 1.5

    assert_equal(c2(psi), expected_c2)
    assert_equal(c3(psi), expected_c3)
