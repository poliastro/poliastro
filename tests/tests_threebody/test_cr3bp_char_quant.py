"""
@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University
"""
from astropy.tests.helper import assert_quantity_allclose

from poliastro.threebody.cr3bp_char_quant import SystemChars


def test_mu(p1,p2,expected_mu):
    "Test cr3bp_char_quant -> SystemChars.mu with expected mu, only compares the value"
    Systemp1p2 = SystemChars(p1,p2)
    
    assert_quantity_allclose(Systemp1p2.mu.value,expected_mu.value,1e-5)
    
def test_lstar(p1,p2,expected_lstar):
    "Test cr3bp_char_quant -> SystemChars.lstar with expected lstar"
    Systemp1p2 = SystemChars(p1,p2)
    
    assert_quantity_allclose(Systemp1p2.lstar,expected_lstar,1e-5)
    
    
def test_tstar(p1,p2,expected_tstar):
    "Test cr3bp_char_quant -> SystemChars.tstar with expected tstar"
    Systemp1p2 = SystemChars(p1,p2)
    
    assert_quantity_allclose(Systemp1p2.tstar,expected_tstar,1e-5)
    
    

