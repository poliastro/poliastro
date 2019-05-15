import numpy as np
import pytest

from poliastro.core import util


def test_rotation_matrix_x():
    result = util.rotation_matrix(0.218, 0)
    expected = np.array([[1.0, 0.0, 0.0], [0.0, 0.97633196, -0.21627739], [0.0, 0.21627739, 0.97633196]])
    assert np.allclose(expected, result)


def test_rotation_matrix_y():
    result = util.rotation_matrix(0.218, 1)
    expected = np.array([[0.97633196, 0.0, 0.21627739], [0.0, 1.0, 0.0], [0.21627739, 0.0, 0.97633196]])
    assert np.allclose(expected, result)


def test_rotation_matrix_z():
    result = util.rotation_matrix(0.218, 2)
    expected = np.array([[0.97633196, -0.21627739, 0.0], [0.21627739, 0.97633196, 0.0], [0.0, 0.0, 1.0]])
    assert np.allclose(expected, result)
