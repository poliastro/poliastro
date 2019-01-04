import pytest

from poliastro.bodies import Jupiter
from poliastro.plotting import OrbitPlotter2D, OrbitPlotter3D
from poliastro.plotting.misc import plot_solar_system
from poliastro.twobody import Orbit


@pytest.mark.parametrize("outer,expected", [(True, 8), (False, 4)])
def test_plot_solar_system_has_expected_number_of_orbits(outer, expected):
    assert len(plot_solar_system(outer).orbits) == expected


@pytest.mark.parametrize(
    "use_3d, plotter_class", [(True, OrbitPlotter3D), (False, OrbitPlotter2D)]
)
def test_plot_solar_system_uses_expected_orbitplotter(use_3d, plotter_class):
    assert isinstance(plot_solar_system(use_3d=use_3d), plotter_class)


def test_redraw_makes_attractor_none():
    op = plot_solar_system()
    op._redraw()
    assert op._attractor_radius is not None


def test_set_frame_plots_same_colors():
    op = plot_solar_system()
    jupiter = Orbit.from_body_ephem(Jupiter)
    op.plot(jupiter)
    colors1 = [orb[2] for orb in op._orbits]
    op.set_frame(*jupiter.pqw())
    colors2 = [orb[2] for orb in op._orbits]
    assert colors1 == colors2
