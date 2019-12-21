from unittest import mock

import pytest
from astropy import units as u

from poliastro.bodies import Earth, Mars, Sun
from poliastro.examples import iss
from poliastro.plotting import OrbitPlotter2D, OrbitPlotter3D


@pytest.mark.parametrize("plotter_class", [OrbitPlotter2D, OrbitPlotter3D])
def test_get_figure_has_expected_properties(plotter_class):
    frame = plotter_class()
    figure = frame.show()

    assert figure.data == ()
    assert figure.layout.autosize is True
    assert "xaxis" in figure.layout
    assert "yaxis" in figure.layout


def test_get_3d_figure_has_expected_properties():
    frame = OrbitPlotter3D()
    figure = frame.show()

    assert figure.data == ()
    assert figure.layout.autosize is True
    assert "xaxis" in figure.layout.scene
    assert "yaxis" in figure.layout.scene
    assert "zaxis" in figure.layout.scene
    assert "aspectmode" in figure.layout.scene


@pytest.mark.parametrize("plotter_class", [OrbitPlotter2D, OrbitPlotter3D])
def test_set_different_attractor_raises_error(plotter_class):
    body1 = Earth

    body2 = Mars

    frame = plotter_class()
    frame.set_attractor(body1)

    with pytest.raises(NotImplementedError) as excinfo:
        frame.set_attractor(body2)
    assert "Attractor has already been set to Earth." in excinfo.exconly()


@pytest.mark.parametrize("plotter_class", [OrbitPlotter2D, OrbitPlotter3D])
def test_plot_sets_attractor(plotter_class):
    frame = plotter_class()
    assert frame._attractor is None

    frame.plot(iss)
    assert frame._attractor == iss.attractor


@pytest.mark.parametrize("plotter_class", [OrbitPlotter2D, OrbitPlotter3D])
def test_plot_appends_data(plotter_class):
    frame = plotter_class()
    assert len(frame.trajectories) == 0

    frame.plot(iss)
    assert len(frame.trajectories) == 1


@pytest.mark.parametrize("plotter_class", [OrbitPlotter2D, OrbitPlotter3D])
def test_plot_trajectory_without_attractor_raises_error(plotter_class):
    frame = plotter_class()

    with pytest.raises(ValueError) as excinfo:
        frame.plot_trajectory({})
    assert (
        "An attractor must be set up first, please use "
        "set_attractor(Major_Body) or plot(orbit)." in excinfo.exconly()
    )


def test_plot_2d_trajectory_without_frame_raises_error():
    frame = OrbitPlotter2D()

    with pytest.raises(ValueError) as excinfo:
        frame.set_attractor(Sun)
        frame.plot_trajectory({})
    assert (
        "A frame must be set up first, please use "
        "set_frame(*orbit.pqw()) or plot(orbit)." in excinfo.exconly()
    )


def test_plot_3d_trajectory_plots_a_trajectory():
    frame = OrbitPlotter3D()
    assert len(frame.trajectories) == 0

    trajectory = Earth.get_mean_orbit().sample()
    frame.set_attractor(Sun)
    frame.plot_trajectory(trajectory)

    assert len(frame.trajectories) == 1
    assert frame._attractor == Sun


def test_plot_2d_trajectory_plots_a_trajectory():
    frame = OrbitPlotter2D()
    assert len(frame.trajectories) == 0

    earth = Earth.get_mean_orbit()
    trajectory = earth.sample()
    frame.set_attractor(Sun)
    frame.set_frame(*earth.pqw())
    frame.plot_trajectory(trajectory)

    assert len(frame.trajectories) == 1
    assert frame._attractor == Sun


@pytest.mark.parametrize("plotter_class", [OrbitPlotter2D, OrbitPlotter3D])
def test_show_calls_prepare_plot(plotter_class):
    with mock.patch.object(plotter_class, "_prepare_plot") as mock_prepare_plot:
        m = plotter_class()
        m.plot(orbit=Earth.get_mean_orbit(), label="Object")
        m.show()

        mock_prepare_plot.assert_called_with()


def test_set_view():
    frame = OrbitPlotter3D()
    frame.set_view(0 * u.deg, 0 * u.deg, 1000 * u.m)
    figure = frame.show()

    eye = figure["layout"]["scene"]["camera"]["eye"]
    assert eye["x"] == 1
    assert eye["y"] == 0
    assert eye["z"] == 0


def test_dark_theme():
    frame = OrbitPlotter3D(dark=True)
    assert frame._layout.template.layout.plot_bgcolor == "rgb(17,17,17)"
