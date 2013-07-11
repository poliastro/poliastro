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

    def test_vectorize_vec(self):
        vec = np.random.rand(3, 10)
        ax = np.array([1, 0, 0])
        angle = np.random.rand()
        rot_vec = rotate(vec, ax, angle)
        assert rot_vec.shape == vec.shape

    def test_vectorize_ang(self):
        vec = np.random.rand(3)
        ax = np.array([1, 0, 0])
        angle = np.random.rand(10)
        rot_vec = rotate(vec, ax, angle)
        assert rot_vec.shape == vec.shape + angle.shape

    def test_data(self):
        N = 10
        idx = 1
        angles = radians(np.arange(N) * 10)
        angle = angles[idx]
        vv = np.array([
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [-1, 0, 0],
            [0, -1, 0],
            [0, 0, -1],
            [1, 1, 0],
            [1, 0, 1],
            [0, 1, 1],
            [-1, 0, -1]
            ]).T
        v = vv[:, idx]
        # vv[idx] rotated angles[idx]
        v0 = np.array([-0.1736,  0.9848,  0.    ])
        assert_array_almost_equal(rotate(vv, 3, angles)[:, 1], v0, decimal=4)
        assert_array_almost_equal(rotate(vv, 3, angle)[:, 1], v0, decimal=4)
        assert_array_almost_equal(rotate(v, 3, angles)[:, 1], v0, decimal=4)
        assert_array_almost_equal(rotate(v, 3, angle), v0, decimal=4)

        assert_array_almost_equal(rotate(vv, 3, angles)[0][1], v0[0], decimal=4)
        assert_array_almost_equal(rotate(vv, 3, angle)[0][1], v0[0], decimal=4)
        assert_array_almost_equal(rotate(v, 3, angles)[0][1], v0[0], decimal=4)
        assert_array_almost_equal(rotate(v, 3, angle)[0], v0[0], decimal=4)


if __name__ == '__main__':
    run_module_suite()
