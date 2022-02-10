import pytest
from astropy import units as u
from matplotlib import pyplot as plt

from poliastro.bodies import Earth
from poliastro.examples import iss
from poliastro.plotting.gabbard import GabbardPlotter
from poliastro.twobody import Orbit


def test_axes_labels_and_title():
    ax = plt.gca()
    gbp = GabbardPlotter(ax=ax)
    ss = iss
    gbp.plot_orbits([ss])

    assert ax.get_xlabel() == "Period (min)"
    assert ax.get_ylabel() == "Altitude (km)"


def test_legend():
    ax = plt.gca()
    gbp = GabbardPlotter(ax=ax)
    ss = iss
    gbp.plot_orbits([ss], label="ISS")
    legend = plt.gca().get_legend()

    ss.epoch.out_subfmt = "date_hm"
    label = f"{ss.epoch.iso} (ISS)"

    assert legend.get_title().get_text() == label


@pytest.mark.mpl_image_compare
def test_static_gabbard_plotting():
    fig, ax = plt.subplots()
    plotter = GabbardPlotter(ax=ax)
    plotter.plot_orbits(debris_orbits_list, label="COSMOS 1408 DEB")

    return fig


@pytest.mark.parametrize(
    "dark, expected_color", [(True, (0.0, 0.0, 0.0, 1.0)), (False, (1.0, 1.0, 1.0, 1))]
)
def test_dark_mode_plots_dark_plot(dark, expected_color):
    gbp = GabbardPlotter(dark=dark)
    assert gbp._ax.get_facecolor() == expected_color


debris_orbits_list = []

debris_orbits_list.append(
    Orbit.from_classical(
        attractor=Earth,
        a=6828.193338988509 * u.km,
        ecc=0.0062140354170737815 * u.one,
        inc=82.69387440482602 * u.deg,
        raan=37.33894561668519 * u.deg,
        argp=200.62393574484153 * u.deg,
        nu=-117.55203086408737 * u.deg,
    )
)

debris_orbits_list.append(
    Orbit.from_classical(
        attractor=Earth,
        a=6821.922877133498 * u.km,
        ecc=0.0023241234515223646 * u.one,
        inc=82.65684766470754 * u.deg,
        raan=36.3401421924121 * u.deg,
        argp=125.29597430617513 * u.deg,
        nu=-151.64963315597913 * u.deg,
    )
)

debris_orbits_list.append(
    Orbit.from_classical(
        attractor=Earth,
        a=6836.825166360441 * u.km,
        ecc=0.004635589373624103 * u.one,
        inc=82.69764910622918 * u.deg,
        raan=36.757861621556614 * u.deg,
        argp=44.219092511353594 * u.deg,
        nu=133.63349740950568 * u.deg,
    )
)

debris_orbits_list.append(
    Orbit.from_classical(
        attractor=Earth,
        a=7117.414488497581 * u.km,
        ecc=0.04540691741651712 * u.one,
        inc=83.07451144780156 * u.deg,
        raan=52.87995597314799 * u.deg,
        argp=190.63045916106168 * u.deg,
        nu=41.306044841636634 * u.deg,
    )
)

debris_orbits_list.append(
    Orbit.from_classical(
        attractor=Earth,
        a=6917.149445506758 * u.km,
        ecc=0.014811336719791061 * u.one,
        inc=82.6218902660939 * u.deg,
        raan=39.58175296436131 * u.deg,
        argp=106.71561062464224 * u.deg,
        nu=-62.66454424955413 * u.deg,
    )
)

debris_orbits_list.append(
    Orbit.from_classical(
        attractor=Earth,
        a=6910.923895002164 * u.km,
        ecc=0.016808366674815105 * u.one,
        inc=82.35369942440258 * u.deg,
        raan=35.60505049154483 * u.deg,
        argp=184.42211913686066 * u.deg,
        nu=45.95875318418421 * u.deg,
    )
)

debris_orbits_list.append(
    Orbit.from_classical(
        attractor=Earth,
        a=6677.880585440615 * u.km,
        ecc=0.0014802577911015675 * u.one,
        inc=82.17121703030627 * u.deg,
        raan=25.699484007134643 * u.deg,
        argp=204.23415215165576 * u.deg,
        nu=-24.40253410856961 * u.deg,
    )
)

debris_orbits_list.append(
    Orbit.from_classical(
        attractor=Earth,
        a=6961.126175813936 * u.km,
        ecc=0.025433234625720075 * u.one,
        inc=82.56625734212793 * u.deg,
        raan=40.073969157354824 * u.deg,
        argp=188.205744852877 * u.deg,
        nu=152.76011672297756 * u.deg,
    )
)

debris_orbits_list.append(
    Orbit.from_classical(
        attractor=Earth,
        a=6883.993995695628 * u.km,
        ecc=0.008022496129527943 * u.one,
        inc=83.2763699331857 * u.deg,
        raan=47.55215625476586 * u.deg,
        argp=114.15854367766484 * u.deg,
        nu=-31.793479924939778 * u.deg,
    )
)

debris_orbits_list.append(
    Orbit.from_classical(
        attractor=Earth,
        a=6798.16535966413 * u.km,
        ecc=0.006497072501655168 * u.one,
        inc=82.61152388022165 * u.deg,
        raan=34.66534775192231 * u.deg,
        argp=349.72219499585407 * u.deg,
        nu=90.16499314429353 * u.deg,
    )
)

debris_orbits_list.append(
    Orbit.from_classical(
        attractor=Earth,
        a=6854.939735394475 * u.km,
        ecc=0.006324699282668879 * u.one,
        inc=82.67294675705102 * u.deg,
        raan=37.411303678153935 * u.deg,
        argp=69.71516857133007 * u.deg,
        nu=-29.51423158409656 * u.deg,
    )
)
debris_orbits_list.append(
    Orbit.from_classical(
        attractor=Earth,
        a=6884.147529089194 * u.km,
        ecc=0.01552728578907267 * u.one,
        inc=82.4114164903979 * u.deg,
        raan=35.234082427651664 * u.deg,
        argp=186.60344739193755 * u.deg,
        nu=-146.9445890288543 * u.deg,
    )
)
debris_orbits_list.append(
    Orbit.from_classical(
        attractor=Earth,
        a=6846.946980540078 * u.km,
        ecc=0.0029371405696242913 * u.one,
        inc=82.6314212152875 * u.deg,
        raan=36.88448947562918 * u.deg,
        argp=134.53438085198738 * u.deg,
        nu=-69.91109773157386 * u.deg,
    )
)
debris_orbits_list.append(
    Orbit.from_classical(
        attractor=Earth,
        a=6914.901161035591 * u.km,
        ecc=0.010664583104155329 * u.one,
        inc=82.31703307692484 * u.deg,
        raan=34.654644407052835 * u.deg,
        argp=144.83292617925179 * u.deg,
        nu=-133.54025144695484 * u.deg,
    )
)
debris_orbits_list.append(
    Orbit.from_classical(
        attractor=Earth,
        a=7040.971575280624 * u.km,
        ecc=0.0333018067425175 * u.one,
        inc=82.50417227979605 * u.deg,
        raan=44.11015739081946 * u.deg,
        argp=133.5425169343891 * u.deg,
        nu=-42.74160359135228 * u.deg,
    )
)
"""Orbit List of COSMOS 1408 Debris Orbits Example

COSMOS 1408 Debris Data Taken from https://celestrak.com/NORAD/

"""
