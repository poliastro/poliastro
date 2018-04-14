import pytest

from unittest import mock

from astropy import units as u

from poliastro.examples import iss
from poliastro.plotting import OrbitPlotter3D


def test_get_figure_has_expected_properties():
    frame = OrbitPlotter3D()
    figure = frame.figure

    assert figure["data"] == [{}]
    assert figure["layout"]["autosize"] is True
    assert "xaxis" in figure["layout"]["scene"]
    assert "yaxis" in figure["layout"]["scene"]
    assert "zaxis" in figure["layout"]["scene"]
    assert "aspectmode" in figure["layout"]["scene"]


def test_set_different_attractor_raises_error():
    body1 = mock.MagicMock()
    body1.name = "1"

    body2 = mock.MagicMock()
    body2.name = "2"

    frame = OrbitPlotter3D()
    frame.set_attractor(body1)

    with pytest.raises(NotImplementedError) as excinfo:
        frame.set_attractor(body2)
    assert "Attractor has already been set to 1." in excinfo.exconly()


def test_plot_sets_attractor():
    frame = OrbitPlotter3D()
    assert frame._attractor is None
    assert frame._attractor_data == {}

    frame.plot(iss)
    assert frame._attractor == iss.attractor
    assert frame._attractor_data["name"] == iss.attractor.name


def test_plot_appends_data():
    frame = OrbitPlotter3D()
    assert len(frame._data) == 0

    frame.plot(iss)
    assert len(frame._data) == 1 + 1


def test_plot_trajectory_without_attractor_raises_error():
    frame = OrbitPlotter3D()

    with pytest.raises(ValueError) as excinfo:
        frame.plot_trajectory({})
    assert ("An attractor must be set up first, please use "
            "set_attractor(Major_Body)." in excinfo.exconly())


def test_set_view():
    frame = OrbitPlotter3D()
    frame.set_view(0 * u.deg, 0 * u.deg, 1000 * u.m)

    eye = frame.figure["layout"]["scene"]["camera"]["eye"]
    assert eye["x"] == 1
    assert eye["y"] == 0
    assert eye["z"] == 0


def test_plot_trajectory_plots_a_trajectory():
    frame = OrbitPlotter3D()
    assert len(frame._data) == 0

    earth = Orbit.from_body_ephem(Earth)
    _, trajectory = earth.sample()
    frame.set_attractor(Sun)
    frame.plot_trajectory(trajectory)
    assert len(frame._data) == 1
    assert frame._attractor == Sun


def test_show_calls_prepare_plot():
    patcher = mock.patch.object(OrbitPlotter3D, '_prepare_plot')
    patched = patcher.start()

    m = OrbitPlotter3D()
    earth = Orbit.from_body_ephem(Earth)
    m.plot(orbit=earth, label="Object")
    m.show()

    assert patched.call_count == 1
    patched.assert_called_with()


def test_savefig_calls_prepare_plot():
    patcher = mock.patch.object(OrbitPlotter3D, '_prepare_plot')
    patched = patcher.start()

    m = OrbitPlotter3D()
    earth = Orbit.from_body_ephem(Earth)
    m.plot(orbit=earth, label="Object")
    m.savefig(filename="a.jpeg")

    assert patched.call_count == 1
    patched.assert_called_with()
