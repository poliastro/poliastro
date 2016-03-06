import pytest

import numpy as np
from numpy.testing import assert_almost_equal
from scipy import special

from poliastro.hyper import hyp2f1b as hyp2f1


@pytest.mark.parametrize("x", np.linspace(0, 1, num=10))
def test_hyp2f1_battin_scalar(x):
    expected_res = special.hyp2f1(3, 1, 5/2, x)

    res = hyp2f1(x)
    assert_almost_equal(res, expected_res)
