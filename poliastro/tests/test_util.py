import numpy as np
from numpy import radians
from numpy.testing import TestCase, assert_array_almost_equal, \
    assert_almost_equal, run_module_suite

from poliastro.util import rotate


class TestRotate(TestCase):
    def test_simple(self):
        vec = np.array([1, 0, 0])
        ax = np.array([0, 0, 1])
        angle = np.pi / 2
        rot_vec = rotate(vec, ax, angle)
        assert_array_almost_equal(rot_vec, np.array([0, 1, 0]))

    def test_vectorize(self):
        vec = np.random.rand(3, 10)
        ax = np.array([1, 0, 0])
        angle = np.random.rand()
        rot_vec = rotate(vec, ax, angle)
        assert rot_vec.shape == vec.shape


if __name__ == '__main__':
    run_module_suite()
