import pytest

from poliastro.plotting import OrbitPlotter2D, OrbitPlotter3D
from poliastro.plotting.misc import plot_solar_system


@pytest.mark.parametrize("outer,expected", [(True, 8), (False, 4)])
def test_plot_solar_system_has_expected_number_of_orbits(outer, expected):
    assert len(plot_solar_system(outer).trajectories) == expected


@pytest.mark.parametrize(
    "use_3d, plotter_class", [(True, OrbitPlotter3D), (False, OrbitPlotter2D)]
)
def test_plot_solar_system_uses_expected_orbitplotter(use_3d, plotter_class):
    assert isinstance(plot_solar_system(use_3d=use_3d, interactive=True), plotter_class)
