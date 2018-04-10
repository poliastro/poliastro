import pytest

from astropy import units as u

from poliastro.examples import iss
from poliastro.plotting import OrbitPlotter2D
from poliastro.bodies import Earth, Mars


def test_get_figure_has_expected_properties():
    frame = OrbitPlotter2D()
    figure = frame.figure

    assert figure["data"] == [{}]
    assert figure["layout"]["autosize"] is True
    assert "xaxis" in figure["layout"]["scene"]
    assert "yaxis" in figure["layout"]["scene"]
    assert "aspectmode" in figure["layout"]["scene"]


def test_set_different_attractor_raises_error():
    body1 = Earth

    body2 = Mars

    frame = OrbitPlotter2D()
    frame.set_attractor(body1)

    with pytest.raises(NotImplementedError) as excinfo:
        frame.set_attractor(body2)
    assert "Attractor has already been set to Earth." in excinfo.exconly()


def test_plot_sets_attractor():
    frame = OrbitPlotter2D()
    assert frame._attractor is None
    assert frame._attractor_data == {}

    frame.plot(iss)
    assert frame._attractor == iss.attractor
    assert frame._attractor_data["name"] == iss.attractor.name


def test_plot_appends_data():
    frame = OrbitPlotter2D()
    assert len(frame._data) == 0

    frame.plot(iss)
    assert len(frame._data) == 1 + 1


def test_plot_trajectory_without_attractor_raises_error():
    frame = OrbitPlotter2D()

    with pytest.raises(ValueError) as excinfo:
        frame.plot_trajectory({})
    assert ("An attractor must be set up first, please use "
            "set_attractor(Major_Body)." in excinfo.exconly())


def test_set_view():
    frame = OrbitPlotter2D()
    frame.set_view(0 * u.deg, 0 * u.deg, 1000 * u.m)

    eye = frame.figure["layout"]["scene"]["camera"]["eye"]
    assert eye["x"] == 1
    assert eye["y"] == 0
