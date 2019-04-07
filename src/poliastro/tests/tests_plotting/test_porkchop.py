from poliastro.plotting.porkchop import porkchop
from poliastro.bodies import Earth, Mars
from poliastro.util import time_range

from astropy import units as u
import pytest
import matplotlib.pyplot as plt


@pytest.mark.mpl_image_compare
def test_porkchop_plotting():

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    launch_span = time_range("2005-04-30", end="2005-10-07")
    arrival_span = time_range("2005-11-16", end="2006-12-21")
    dv_dpt, dv_arr, c3dpt, c3arr, tof = porkchop(Earth, Mars, launch_span, arrival_span, ax=ax)

    return fig

