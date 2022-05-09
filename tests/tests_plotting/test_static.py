import sys

import pytest
from astropy import units as u
from astropy.coordinates import CartesianDifferential, CartesianRepresentation
from astropy.time import Time
from matplotlib import pyplot as plt

from poliastro.bodies import Earth, Jupiter, Mars, Sun
from poliastro.constants import J2000_TDB
from poliastro.ephem import Ephem
from poliastro.examples import churi, iss, molniya
from poliastro.frames import Planes
from poliastro.maneuver import Maneuver
from poliastro.plotting.static import StaticOrbitPlotter
from poliastro.twobody import Orbit
from poliastro.util import time_range


def test_axes_labels_and_title():
    ax = plt.gca()
    op = StaticOrbitPlotter(ax)
    ss = iss
    op.plot(ss)

    assert ax.get_xlabel() == "$x$ (km)"
    assert ax.get_ylabel() == "$y$ (km)"


def test_number_of_lines_for_osculating_orbit():
    op1 = StaticOrbitPlotter()
    ss = iss

    l1 = op1.plot(ss)

    assert len(l1) == 2


def test_legend():
    op = StaticOrbitPlotter()
    ss = iss
    op.plot(ss, label="ISS")
    legend = plt.gca().get_legend()

    ss.epoch.out_subfmt = "date_hm"
    label = f"{ss.epoch.iso} (ISS)"

    assert legend.get_texts()[0].get_text() == label


def test_color():
    op = StaticOrbitPlotter()
    ss = iss
    c = "#FF0000"
    op.plot(ss, label="ISS", color=c)
    ax = plt.gca()

    assert ax.get_legend().get_lines()[0].get_c() == c
    for element in ax.get_lines():
        assert element.get_c() == c


def test_plot_trajectory_sets_label():
    expected_label = "67P"

    op = StaticOrbitPlotter()
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
    op = StaticOrbitPlotter(dark=dark)
    assert op._ax.get_facecolor() == expected_color


def test_redraw_makes_attractor_none():
    # TODO: Review
    op = StaticOrbitPlotter()
    op._redraw()
    assert op._attractor_radius is not None


def test_set_frame_plots_same_colors():
    # TODO: Review
    op = StaticOrbitPlotter()
    op.plot_body_orbit(Jupiter, J2000_TDB)
    colors1 = [orb[2] for orb in op.trajectories]
    op.set_body_frame(Jupiter)
    colors2 = [orb[2] for orb in op.trajectories]
    assert colors1 == colors2


def test_redraw_keeps_trajectories():
    # See https://github.com/poliastro/poliastro/issues/518
    op = StaticOrbitPlotter()
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

    op = StaticOrbitPlotter(plane=Planes.EARTH_ECLIPTIC)
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
    plotter = StaticOrbitPlotter(ax=ax)
    plotter.plot(iss)

    return fig


@pytest.mark.mpl_image_compare
def test_basic_trajectory_plotting():
    fig, ax = plt.subplots()
    plotter = StaticOrbitPlotter(ax=ax)
    plotter.set_attractor(Earth)
    plotter.set_orbit_frame(iss)
    plotter.plot_trajectory(iss.sample())

    return fig


@pytest.mark.mpl_image_compare
def test_basic_orbit_and_trajectory_plotting():
    fig, ax = plt.subplots()
    plotter = StaticOrbitPlotter(ax=ax)
    plotter.plot(iss)
    plotter.plot_trajectory(molniya.sample(), label="Molniya")

    return fig


@pytest.mark.mpl_image_compare
def test_trail_plotting():
    fig, ax = plt.subplots()
    plotter = StaticOrbitPlotter(ax=ax)
    plotter.plot(iss, trail=True)

    return fig


@pytest.mark.mpl_image_compare
def test_plot_different_planes():
    fig, ax = plt.subplots()
    plotter = StaticOrbitPlotter(ax=ax)
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
    plotter = StaticOrbitPlotter(ax=ax)
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
    plotter = StaticOrbitPlotter(ax=ax)
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

    op = StaticOrbitPlotter()
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
    ss_i = Orbit.from_classical(
        attractor=Earth,
        a=Earth.R + alt_i,
        ecc=0 * u.one,
        inc=_a,
        raan=_a,
        argp=_a,
        nu=_a,
    )

    # Create the maneuver
    man = Maneuver.hohmann(ss_i, Earth.R + alt_f)

    # Plot the maneuver
    fig, ax = plt.subplots()
    plotter = StaticOrbitPlotter(ax=ax)
    plotter.plot(ss_i, label="Initial orbit", color="blue")
    plotter.plot_maneuver(ss_i, man, label="Hohmann maneuver", color="red")

    return fig
