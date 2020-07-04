import pytest
from astropy import units as u
from matplotlib import pyplot as plt

from poliastro.bodies import Earth, Mars, Venus
from poliastro.plotting._base import BODY_COLORS
from poliastro.plotting.tisserand import TisserandKind, TisserandPlotter


@pytest.mark.mpl_image_compare
def test_tisserand_plotting():
    fig, ax = plt.subplots()

    # Build custom axis
    fig, ax = plt.subplots(1, 1, figsize=(15, 7))
    ax.set_title("Energy Tisserand for Venus, Earth and Mars")
    ax.set_xlabel("$R_{p} [AU]$")
    ax.set_ylabel("Heliocentric Energy [km2 / s2]")
    ax.set_xscale("log")
    ax.set_xlim(10 ** -0.4, 10 ** 0.15)
    ax.set_ylim(-700, 0)

    # Generate a Tisserand plotter
    tp = TisserandPlotter(axes=ax, kind=TisserandKind.ENERGY)

    # Plot Tisserand lines within 1km/s and 10km/s
    for planet in [Venus, Earth, Mars]:
        ax = tp.plot(planet, (1, 14) * u.km / u.s, num_contours=14)

    # Let us label previous figure
    tp.ax.text(0.70, -650, "Venus", color=BODY_COLORS["Venus"])
    tp.ax.text(0.95, -500, "Earth", color=BODY_COLORS["Earth"])
    tp.ax.text(1.35, -350, "Mars", color=BODY_COLORS["Mars"])

    return fig
