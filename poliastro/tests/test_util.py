import numpy as np
from numpy import radians
from numpy.testing import TestCase, assert_array_almost_equal, \
    assert_almost_equal, run_module_suite

from poliastro.util import rotate, direct_angles


class TestRotate(TestCase):
    def test_rotation_around_coordinate_axis(self):
        vec = np.array([1, 0, 0])
        ax = 2
        angle = np.pi / 2
        rot_vec = rotate(vec, ax, angle)
        assert_array_almost_equal(rot_vec, np.array([0, 0, -1]))

    def test_rotation_around_arbitrary_axis(self):
        vec = np.array([1, 0, 0])
        ax = np.array([0, 0, 1])
        angle = np.pi / 2
        rot_vec = rotate(vec, ax, angle)
        assert_array_almost_equal(rot_vec, np.array([0, 1, 0]))

    def test_vectorize_vec(self):
        vec = np.random.rand(3, 10)
        ax = 1
        angle = np.random.rand()
        rot_vec = rotate(vec, ax, angle)
        assert rot_vec.shape == vec.shape


class TestDirectAngles(TestCase):
    def test_data(self):
        rd = radians
        assert direct_angles(0, rd(-1)) == (0, rd(359))
        assert direct_angles(0, rd(-361)) == (0, rd(359))
        assert direct_angles(rd(360), rd(-361)) == (rd(360), rd(719))


if __name__ == '__main__':
    run_module_suite()
