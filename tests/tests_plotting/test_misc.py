import pytest
from matplotlib import pyplot as plt

from poliastro.plotting.misc import plot_solar_system
from poliastro.plotting.orbit.backends import SUPPORTED_ORBIT_PLOTTER_BACKENDS


@pytest.mark.parametrize("outer, expected", [(True, 8), (False, 4)])
def test_plot_solar_system_has_expected_number_of_orbits(outer, expected):
    assert len(plot_solar_system(outer=outer).trajectories) == expected


@pytest.mark.parametrize("backend_name", SUPPORTED_ORBIT_PLOTTER_BACKENDS)
def test_plot_solar_system_uses_expected_orbitplotter(backend_name):
    assert isinstance(
        plot_solar_system(backend_name=backend_name).backend,
        SUPPORTED_ORBIT_PLOTTER_BACKENDS[backend_name],
    )


@pytest.mark.mpl_image_compare
def test_plot_inner_solar_system_using_matplotlib2D_backend(earth_perihelion):
    plot_solar_system(
        epoch=earth_perihelion,
        outer=False,
        backend_name="matplotlib2D",
    )
    return plt.gcf()


@pytest.mark.mpl_image_compare
def test_plot_outer_solar_system_using_matplotlib2D_backend(earth_perihelion):
    plot_solar_system(
        epoch=earth_perihelion,
        outer=True,
        backend_name="matplotlib2D",
    )
    return plt.gcf()
