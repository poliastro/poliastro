import sys

import pytest
from astropy import time, units as u
from astropy.coordinates import CartesianDifferential, CartesianRepresentation
from astropy.time import Time
from matplotlib import pyplot as plt


from poliastro.bodies import Earth, Jupiter, Mars, Sun
from poliastro.constants import J2000_TDB
from poliastro.ephem import Ephem
from poliastro.examples import churi, iss, molniya
from poliastro.frames import Planes
from poliastro.maneuver import Maneuver
from poliastro.twobody import Orbit
from poliastro.plotting import OrbitPlotter
from poliastro.util import time_range


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
    assert "Attractor has already been set to Earth" in excinfo.exconly()


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
        frame._frame = 1  # Set it to something not None to skip frame check
        frame.plot_trajectory({})
    assert (
        "An attractor must be set up first, please use "
        "set_attractor(Major_Body) or plot(orbit)" in excinfo.exconly()
    )


def test_plot_2d_trajectory_without_frame_raises_error():
    frame = OrbitPlotter2D()

    with pytest.raises(ValueError) as excinfo:
        frame.set_attractor(Sun)
        frame.plot_trajectory({})
    assert (
        "A frame must be set up first, please use "
        "set_orbit_frame(orbit) or plot(orbit)" in excinfo.exconly()
    )


def test_ephem_without_frame_raises_error():
    epochs = time.Time("2020-04-29 10:43", scale="tdb")
    earth = Ephem.from_body(Earth, epochs)
    plotter = OrbitPlotter2D()

    with pytest.raises(ValueError) as excinfo:
        plotter.set_attractor(Sun)
        plotter.plot_ephem(earth)
    assert (
        "A frame must be set up first, please use "
        "set_orbit_frame(orbit) or plot(orbit)" in excinfo.exconly()
    )


def test_plot_3d_trajectory_plots_a_trajectory():
    frame = OrbitPlotter3D()
    assert len(frame.trajectories) == 0

    trajectory = churi.sample()
    frame.set_attractor(Sun)
    frame.plot_trajectory(trajectory)

    assert len(frame.trajectories) == 1
    assert frame._attractor == Sun


def test_plot_2d_trajectory_plots_a_trajectory():
    frame = OrbitPlotter2D()
    assert len(frame.trajectories) == 0

    trajectory = churi.sample()
    frame.set_attractor(Sun)
    frame.set_orbit_frame(churi)
    frame.plot_trajectory(trajectory)

    assert len(frame.trajectories) == 1
    assert frame._attractor == Sun


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



def test_axes_labels_and_title():
    ax = plt.gca()
    op = OrbitPlotter(ax)
    ss = iss
    op.plot(ss)

    assert ax.get_xlabel() == "$x$ (km)"
    assert ax.get_ylabel() == "$y$ (km)"


def test_number_of_lines_for_osculating_orbit():
    op1 = OrbitPlotter(backend_name="matplotlib2D")
    ss = iss

    l1 = op1.plot(ss)

    assert len(l1) == 2


def test_legend():
    op = OrbitPlotter(backend_name="matplotlib2D")
    ss = iss
    op.plot(ss, label="ISS")
    legend = plt.gca().get_legend()

    ss.epoch.out_subfmt = "date_hm"
    label = f"{ss.epoch.iso} (ISS)"

    assert legend.get_texts()[0].get_text() == label


def test_color():
    op = OrbitPlotter(backend_name="matplotlib2D")
    ss = iss
    c = "#FF0000"
    op.plot(ss, label="ISS", color=c)
    ax = plt.gca()

    assert ax.get_legend().get_lines()[0].get_c() == c
    for element in ax.get_lines():
        assert element.get_c() == c


def test_plot_trajectory_sets_label():
    expected_label = "67P"

    op = OrbitPlotter(backend_name="matplotlib2D")
    trajectory = churi.sample()
    op.plot_body_orbit(Mars, J2000_TDB, label="Mars")

    op.plot_trajectory(trajectory, label=expected_label)

    legend = plt.gca().get_legend()
    assert legend.get_texts()[1].get_text() == expected_label


@pytest.mark.parametrize(
    "dark, expected_color",
    [(True, (0.0, 0.0, 0.0, 1.0)), (False, (1.0, 1.0, 1.0, 1))],
)
def test_dark_mode_plots_dark_plot(dark, expected_color):
    op = OrbitPlotter(dark=dark)
    assert op._ax.get_facecolor() == expected_color


def test_redraw_makes_attractor_none():
    # TODO: Review
    op = OrbitPlotter(backend_name="matplotlib2D")
    op._redraw()
    assert op._attractor_radius is not None


def test_set_frame_plots_same_colors():
    # TODO: Review
    op = OrbitPlotter(backend_name="matplotlib2D")
    op.plot_body_orbit(Jupiter, J2000_TDB)
    colors1 = [orb[2] for orb in op.trajectories]
    op.set_body_frame(Jupiter)
    colors2 = [orb[2] for orb in op.trajectories]
    assert colors1 == colors2


def test_redraw_keeps_trajectories():
    # See https://github.com/poliastro/poliastro/issues/518
    op = OrbitPlotter(backend_name="matplotlib2D")
    trajectory = churi.sample()
    op.plot_body_orbit(Mars, J2000_TDB, label="Mars")
    op.plot_trajectory(trajectory, label="67P")

    assert len(op.trajectories) == 2

    op.set_body_frame(Mars)

    assert len(op.trajectories) == 2


def test_plot_ephem_different_plane_raises_error():
    unused_epochs = Time.now().reshape(-1)
    unused_coordinates = CartesianRepresentation(
        [(1, 0, 0)] * u.au,
        xyz_axis=1,
        differentials=CartesianDifferential(
            [(0, 1, 0)] * (u.au / u.day), xyz_axis=1
        ),
    )

    op = OrbitPlotter(backend_name="matplotlib2D", plane=Planes.EARTH_ECLIPTIC)
    op.set_attractor(Sun)
    op.set_body_frame(Earth)
    with pytest.raises(ValueError) as excinfo:
        op.plot_ephem(
            Ephem(unused_epochs, unused_coordinates, Planes.EARTH_EQUATOR)
        )

    assert (
        "sample the ephemerides using a different plane or create a new plotter"
        in excinfo.exconly()
    )


@pytest.mark.mpl_image_compare
def test_basic_plotting():
    fig, ax = plt.subplots()
    plotter = OrbitPlotter(scene=ax)
    plotter.plot(iss)

    return fig


@pytest.mark.mpl_image_compare
def test_basic_trajectory_plotting():
    fig, ax = plt.subplots()
    plotter = OrbitPlotter(scene=ax)
    plotter.set_attractor(Earth)
    plotter.set_orbit_frame(iss)
    plotter.plot_trajectory(iss.sample())

    return fig


@pytest.mark.mpl_image_compare
def test_basic_orbit_and_trajectory_plotting():
    fig, ax = plt.subplots()
    plotter = OrbitPlotter(scene=ax)
    plotter.plot(iss)
    plotter.plot_trajectory(molniya.sample(), label="Molniya")

    return fig


@pytest.mark.mpl_image_compare
def test_trail_plotting():
    fig, ax = plt.subplots()
    plotter = OrbitPlotter(scene=ax)
    plotter.plot(iss, trail=True)

    return fig


@pytest.mark.mpl_image_compare
def test_plot_different_planes():
    fig, ax = plt.subplots()
    plotter = OrbitPlotter(scene=ax)
    plotter.plot(iss)
    plotter.plot(molniya.change_plane(Planes.EARTH_ECLIPTIC))

    return fig


@pytest.mark.mpl_image_compare
def test_body_plotting(earth_perihelion):
    Earth.plot(earth_perihelion)

    return plt.gcf()


@pytest.mark.mpl_image_compare
@pytest.mark.remote_data
def test_plot_ephem_epoch():
    epoch = Time("2020-02-14 00:00:00")
    ephem = Ephem.from_horizons(
        "2020 CD3",
        time_range(
            Time("2020-02-13 12:00:00"), end=Time("2020-02-14 12:00:00")
        ),
        attractor=Earth,
    )

    fig, ax = plt.subplots()
    plotter = OrbitPlotter(scene=ax)
    plotter.set_attractor(Earth)
    plotter.set_orbit_frame(Orbit.from_ephem(Earth, ephem, epoch))

    plotter.plot_ephem(ephem, epoch, label="2020 CD3 Minimoon", color="k")

    return fig


@pytest.mark.mpl_image_compare
@pytest.mark.remote_data
def test_plot_ephem_no_epoch():
    epoch = Time("2020-02-14 00:00:00")
    ephem = Ephem.from_horizons(
        "2020 CD3",
        time_range(
            Time("2020-02-13 12:00:00"), end=Time("2020-02-14 12:00:00")
        ),
        attractor=Earth,
    )

    fig, ax = plt.subplots()
    plotter = OrbitPlotter(scene=ax)
    plotter.set_attractor(Earth)
    plotter.set_orbit_frame(Orbit.from_ephem(Earth, ephem, epoch))

    plotter.plot_ephem(ephem, label="2020 CD3 Minimoon", color="k")

    return fig


def test_body_frame_raises_warning_if_time_is_not_tdb_with_proper_time(
    recwarn,
):
    from poliastro.warnings import TimeScaleWarning

    body = Jupiter
    epoch = Time("2017-09-29 07:31:26", scale="utc")
    expected_epoch_string = "2017-09-29 07:32:35.182"  # epoch.tdb.value

    op = OrbitPlotter(backend_name="matplotlib2D")
    op.set_body_frame(body, epoch)

    w = recwarn.pop(TimeScaleWarning)

    assert expected_epoch_string in str(w.message)


@pytest.mark.xfail(
    sys.maxsize < 2**32, reason="not supported for 32 bit systems"
)
@pytest.mark.mpl_image_compare
def test_plot_maneuver():
    # Data from Vallado, example 6.1
    alt_i = 191.34411 * u.km
    alt_f = 35781.34857 * u.km
    _a = 0 * u.deg
    orb_i = Orbit.from_classical(
        attractor=Earth,
        a=Earth.R + alt_i,
        ecc=0 * u.one,
        inc=_a,
        raan=_a,
        argp=_a,
        nu=_a,
    )

    # Create the maneuver
    man = Maneuver.hohmann(orb_i, Earth.R + alt_f)

    # Plot the maneuver
    fig, ax = plt.subplots()
    plotter = OrbitPlotter(scene=ax, backend_name="matplotlib2D")
    plotter.plot(orb_i, label="Initial orbit", color="blue")
    plotter.plot_maneuver(orb_i, man, label="Hohmann maneuver", color="red")

    return fig
