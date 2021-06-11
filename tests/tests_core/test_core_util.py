from functools import partial

import numpy as np
import pytest
from hypothesis import given, settings, strategies as st

from poliastro.core.util import alinspace, rotation_matrix


def test_rotation_matrix_x():
    result = rotation_matrix(0.218, 0)
    expected = np.array(
        [[1.0, 0.0, 0.0], [0.0, 0.97633196, -0.21627739], [0.0, 0.21627739, 0.97633196]]
    )
    assert np.allclose(expected, result)


def test_rotation_matrix_y():
    result = rotation_matrix(0.218, 1)
    expected = np.array(
        [[0.97633196, 0.0, 0.21627739], [0.0, 1.0, 0.0], [0.21627739, 0.0, 0.97633196]]
    )
    assert np.allclose(expected, result)


def test_rotation_matrix_z():
    result = rotation_matrix(0.218, 2)
    expected = np.array(
        [[0.97633196, -0.21627739, 0.0], [0.21627739, 0.97633196, 0.0], [0.0, 0.0, 1.0]]
    )
    assert np.allclose(expected, result)


angles = partial(st.floats, min_value=-2 * np.pi, max_value=2 * np.pi)


@settings(deadline=None)
@given(
    x=angles(),
    y=st.one_of(angles(), st.none()),
)
def test_alinspace_is_always_increasing(x, y):
    result = alinspace(x, y)

    assert (np.diff(result) >= 0).all()


@settings(deadline=None)
@given(
    x=st.floats(min_value=-np.pi, max_value=np.pi),
    y=st.one_of(st.floats(min_value=-np.pi, max_value=np.pi), st.none()),
)
def test_alinspace_is_always_increasing_with_angles_inside_anomaly_range(x, y):
    result = alinspace(x, y)

    assert (np.diff(result) >= 0).all()


@settings(deadline=None)
@given(x=angles())
def test_alinspace_no_max_value_uses_full_circle(x):
    result = alinspace(x)

    assert result.max() - result.min() == pytest.approx(2 * np.pi)
