import pytest

from poliastro.examples import iss
from poliastro.plotting import AnimatedOrbitPlotter


def test_type_of_anim():
    op = AnimatedOrbitPlotter()
    ss = iss
    k = op.anim(ss)
    assert type(k) == tuple


@pytest.mark.parametrize(
    "dark, expected_color", [(True, (0.0, 0.0, 0.0, 1.0)), (False, (1.0, 1.0, 1.0, 1))]
)
def test_dark_mode_plots_dark_plot(dark, expected_color):
    op = AnimatedOrbitPlotter(dark=dark)
    assert op._ax.get_facecolor() == expected_color
