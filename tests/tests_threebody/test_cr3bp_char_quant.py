"""
@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University
"""
import pytest
from astropy.tests.helper import assert_quantity_allclose

from poliastro.bodies import Earth, Moon
from poliastro.core.threebody.cr3bp_quantities_calculations import (
    calculate_mu,
    calculate_tstar,
)
from poliastro.threebody.cr3bp_char_quant import SystemChars


@pytest.mark.parametrize(
    "p1, p2, expected_mu",
    [
        (Earth, Moon, calculate_mu(Earth.k, Moon.k)),
        # Compares value from SystemChars class with the preivously tested expected_mu()
    ],
)
def test_SystemChars_mu(p1, p2, expected_mu):
    "Test cr3bp_char_quant -> SystemChars.mu with expected mu"
    Systemp1p2 = SystemChars.from_primaries(p1, p2)

    assert_quantity_allclose(Systemp1p2.mu, expected_mu, 1e-5)


@pytest.mark.parametrize(
    "p1, p2, expected_lstar",
    [
        (Earth, Moon, Moon.mean_a),
    ],
)
def test_lstar(p1, p2, expected_lstar):
    "Test cr3bp_char_quant -> SystemChars.lstar with expected lstar"
    Systemp1p2 = SystemChars.from_primaries(p1, p2)

    assert_quantity_allclose(Systemp1p2.lstar, expected_lstar, 1e-5)


@pytest.mark.parametrize(
    "p1, p2, expected_tstar",
    [
        (Earth, Moon, calculate_tstar(Earth.k, Moon.k, Moon.mean_a)),
        # Compares value from SystemChars class with the preivously tested expected_tstar()
    ],
)
def test_tstar(p1, p2, expected_tstar):
    "Test cr3bp_char_quant -> SystemChars.tstar with expected tstar"
    Systemp1p2 = SystemChars.from_primaries(p1, p2)

    assert_quantity_allclose(Systemp1p2.tstar, expected_tstar, 1e-5)
