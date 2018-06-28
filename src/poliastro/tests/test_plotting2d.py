import os

import pytest

import tempfile

from astropy import units as u

from unittest import mock

from poliastro.examples import iss
from poliastro.plotting import OrbitPlotter2D
from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody.orbit import Orbit


def test_get_figure_has_expected_properties():
    frame = OrbitPlotter2D()
    figure = frame.figure

    assert figure["data"] == [{}]
    assert figure["layout"]["autosize"] is True
    assert "xaxis" in figure["layout"]
    assert "yaxis" in figure["layout"]


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


def test_plot_trajectory_plots_a_trajectory():
    frame = OrbitPlotter2D()
    assert len(frame._data) == 0

    earth = Orbit.from_body_ephem(Earth)
    _, trajectory = earth.sample()
    frame.set_attractor(Sun)
    frame.plot_trajectory(trajectory)
    assert len(frame._data) == 1
    assert frame._attractor == Sun


def test_show_calls_prepare_plot():
    patcher = mock.patch.object(OrbitPlotter2D, '_prepare_plot')
    patched = patcher.start()

    m = OrbitPlotter2D()
    earth = Orbit.from_body_ephem(Earth)
    m.plot(orbit=earth, label="Obj")
    m.show()

    assert patched.call_count == 1
    patched.assert_called_with()


def test_savefig_calls_prepare_plot():
    patcher = mock.patch.object(OrbitPlotter2D, '_prepare_plot')
    patched = patcher.start()

    m = OrbitPlotter2D()
    earth = Orbit.from_body_ephem(Earth)
    m.plot(orbit=earth, label="Obj")
    with tempfile.NamedTemporaryFile() as fp:
        m.savefig(filename=fp.name+".jpeg")

    assert patched.call_count == 1
    patched.assert_called_with()
