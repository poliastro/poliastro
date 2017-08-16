# coding: utf-8
import pytest

import astropy.units as u

import matplotlib.pyplot as plt

from poliastro.examples import iss

from poliastro.plotting import OrbitPlotter


def test_orbitplotter_has_axes():
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


def test_axes_labels_and_title():
    ax = plt.gca()
    op = OrbitPlotter(ax)
    ss = iss
    op.plot(ss)

    assert ax.get_xlabel() == "$x$ (km)"
    assert ax.get_ylabel() == "$y$ (km)"
    assert ax.get_title() == str(ss.epoch.iso)


def test_number_of_lines_for_osculating_orbit():
    _, (ax1, ax2) = plt.subplots(ncols=2)
    op1 = OrbitPlotter(ax1)
    op2 = OrbitPlotter(ax2)
    ss = iss

    l1 = op1.plot(ss, osculating=False)
    l2 = op2.plot(ss, osculating=True)

    assert len(l1) == 1
    assert len(l2) == 2
