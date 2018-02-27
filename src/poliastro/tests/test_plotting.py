import pytest

import astropy.units as u

import matplotlib.pyplot as plt

from poliastro.examples import iss

from poliastro.plotting import OrbitPlotter, OrbitPlotter3D, plot_solar_system


def test_orbitplotter_has_axes():
    ax = "Unused axes"
    op = OrbitPlotter(ax)
    assert op.ax is ax


def test_set_frame():
    op = OrbitPlotter()
    p = [1, 0, 0] * u.one
    q = [0, 1, 0] * u.one
    w = [0, 0, 1] * u.one
    op.set_frame(p, q, w)

    assert op._frame == (p, q, w)


def test_axes_labels_and_title():
    ax = plt.gca()
    op = OrbitPlotter(ax)
    ss = iss
    op.plot(ss)

    assert ax.get_xlabel() == "$x$ (km)"
    assert ax.get_ylabel() == "$y$ (km)"


def test_number_of_lines_for_osculating_orbit():
    op1 = OrbitPlotter()
    ss = iss

    l1 = op1.plot(ss)

    assert len(l1) == 2


def test_legend():
    op = OrbitPlotter()
    ss = iss
    op.plot(ss, label='ISS')
    legend = plt.gca().get_legend()

    ss.epoch.out_subfmt = 'date_hm'
    label = '{} ({})'.format(ss.epoch.iso, 'ISS')

    assert legend.get_texts()[0].get_text() == label


def test_color():
    op = OrbitPlotter()
    ss = iss
    c = "#FF0000"
    op.plot(ss, label='ISS', color=c)
    ax = plt.gca()

    assert ax.get_legend().get_lines()[0].get_c() == c
    for element in ax.get_lines():
        assert element.get_c() == c


@pytest.mark.parametrize("outer,expected,function", [
    (True, 8, OrbitPlotter),
    (False, 4, OrbitPlotter),
    (False, 4, OrbitPlotter3D),
])
def test_plot_solar_system(outer, expected, function):
    if function == OrbitPlotter:
        assert len(plot_solar_system(outer, function=OrbitPlotter).orbits) == expected
    assert isinstance(plot_solar_system(function=OrbitPlotter), OrbitPlotter)
    assert isinstance(plot_solar_system(function=OrbitPlotter3D), OrbitPlotter3D)
