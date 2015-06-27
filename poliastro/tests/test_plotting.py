# coding: utf-8
import pytest

import astropy.units as u

import matplotlib
matplotlib.use("Agg", warn=False)  # use Agg backend for these tests

from poliastro import plotting
from poliastro.plotting import OrbitPlotter


def test_OrbitPlotter_has_axes():
    ax = "Unused axes"
    op = OrbitPlotter(ax)
    assert op.ax is ax


def test_set_frame_raises_error_if_frame_exists():
    op = OrbitPlotter()
    p = [1, 0, 0] * u.one
    q = [0, 1, 0] * u.one
    w = [0, 0, 1] * u.one
    op.set_frame(p, q, w)
    with pytest.raises(NotImplementedError):
        op.set_frame(p, q, w)
