from matplotlib import pyplot as plt
import pytest

from poliastro.plotting.misc import plot_solar_system
from poliastro.plotting.orbit.backends import DEFAULT_ORBIT_PLOTTER_BACKENDS
from poliastro.plotting.orbit.backends._base import OrbitPlotterBackend


@pytest.mark.parametrize("outer, expected", [(True, 8), (False, 4)])
def test_plot_solar_system_has_expected_number_of_orbits(outer, expected):
    assert len(plot_solar_system(outer=outer).trajectories) == expected


@pytest.mark.parametrize("Backend", DEFAULT_ORBIT_PLOTTER_BACKENDS.values())
def test_plot_solar_system_uses_expected_orbitplotter(Backend):
    assert isinstance(
        plot_solar_system(backend=Backend()).backend,
        OrbitPlotterBackend,
    )


@pytest.mark.mpl_image_compare
def test_plot_inner_solar_system_using_matplotlib2D_backend(earth_perihelion):
    plot_solar_system(
        epoch=earth_perihelion,
        outer=False,
    )
    return plt.gcf()


@pytest.mark.mpl_image_compare
def test_plot_outer_solar_system_using_matplotlib2D_backend(earth_perihelion):
    plot_solar_system(
        epoch=earth_perihelion,
        outer=True,
    )
    return plt.gcf()
