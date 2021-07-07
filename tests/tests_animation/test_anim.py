import matplotlib
import pytest

from poliastro.examples import churi, iss
from poliastro.plotting import AnimatedOrbitPlotter


def test_type_of_anim():
    op = AnimatedOrbitPlotter()
    ss = iss
    k = op.anim(ss)
    assert type(k) == matplotlib.animation.FuncAnimation


def test_type_of_anim_larger_orbits():
    op = AnimatedOrbitPlotter()
    ss = churi
    k = op.anim(ss)
    assert type(k) == matplotlib.animation.FuncAnimation


@pytest.mark.parametrize(
    "dark, expected_color", [(True, (0.0, 0.0, 0.0, 1.0)), (False, (1.0, 1.0, 1.0, 1))]
)
def test_dark_mode_plots_dark_plot(dark, expected_color):
    op = AnimatedOrbitPlotter(dark=dark)
    assert op._ax.get_facecolor() == expected_color
