import pytest

import tempfile

from unittest import mock

from poliastro.examples import iss
from poliastro.plotting import OrbitPlotter2D
from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody.orbit import Orbit

@pytest.fixture(scope="module")
def op():
    return OrbitPlotter2D()

def test_get_figure_has_expected_properties(op):
    figure = op.figure

    assert figure["data"] == [{}]
    assert figure["layout"]["autosize"] is True
    assert "xaxis" in figure["layout"]
    assert "yaxis" in figure["layout"]


def test_set_different_attractor_raises_error(op):
    body1 = Earth

    body2 = Mars

    op.set_attractor(body1)

    with pytest.raises(NotImplementedError) as excinfo:
        op.set_attractor(body2)
    assert "Attractor has already been set to Earth." in excinfo.exconly()


def test_plot_sets_attractor(op):
    assert op._attractor is None
    assert op._attractor_data == {}

    op.plot(iss)
    assert op._attractor == iss.attractor
    assert op._attractor_data["name"] == iss.attractor.name


def test_plot_appends_data(op):
    assert len(op._data) == 0

    op.plot(iss)
    assert len(op._data) == 1 + 1


def test_plot_trajectory_without_attractor_raises_error(op):
    with pytest.raises(ValueError) as excinfo:
        op.plot_trajectory({})
    assert ("An attractor must be set up first, please use "
            "set_attractor(Major_Body)." in excinfo.exconly())


def test_plot_trajectory_plots_a_trajectory(op):
    assert len(op._data) == 0

    earth = Orbit.from_body_ephem(Earth)
    trajectory = earth.sample()
    op.set_attractor(Sun)
    op.plot_trajectory(trajectory)
    assert len(op._data) == 1
    assert op._attractor == Sun


@mock.patch("poliastro.plotting.iplot")
@mock.patch.object(OrbitPlotter2D, '_prepare_plot')
def test_show_calls_prepare_plot(op, mock_prepare_plot, mock_iplot):
    earth = Orbit.from_body_ephem(Earth)
    op.plot(orbit=earth, label="Obj")
    op.show()

    assert mock_iplot.call_count == 1
    mock_prepare_plot.assert_called_once_with()


@mock.patch("poliastro.plotting.export")
@mock.patch.object(OrbitPlotter2D, '_prepare_plot')
def test_savefig_calls_prepare_plot(op, mock_prepare_plot, mock_export):
    earth = Orbit.from_body_ephem(Earth)
    op.plot(orbit=earth, label="Obj")
    with tempfile.NamedTemporaryFile() as fp:
        m.savefig(filename=fp.name + ".jpeg")

    assert mock_export.call_count == 1
    mock_prepare_plot.assert_called_once_with()
