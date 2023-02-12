from matplotlib import pyplot as plt
import pytest

from poliastro.bodies import Earth, Mars
from poliastro.plotting.porkchop import PorkchopPlotter
from poliastro.util import time_range


@pytest.mark.mpl_image_compare
def test_porkchop_plotting():
    fig, ax = plt.subplots()

    launch_span = time_range("2005-04-30", end="2005-10-07")
    arrival_span = time_range("2005-11-16", end="2006-12-21")
    porkchop_plot = PorkchopPlotter(
        Earth, Mars, launch_span, arrival_span, ax=ax
    )
    dv_dpt, dv_arr, c3dpt, c3arr, tof = porkchop_plot.porkchop()

    return fig
