import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy import special

from poliastro._math.special import hyp2f1b as hyp2f1


@pytest.mark.parametrize("x", np.linspace(0, 1, num=11))
def test_hyp2f1_battin_scalar(x):
    expected_res = special.hyp2f1(3, 1, 5 / 2, x)

    res = hyp2f1(x)
    assert_allclose(res, expected_res)
