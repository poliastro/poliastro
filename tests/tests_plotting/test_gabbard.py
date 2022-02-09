import pytest
from matplotlib import pyplot as plt

from poliastro.examples import debris_orbits_list, iss
from poliastro.plotting.gabbard import GabbardPlotter


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
