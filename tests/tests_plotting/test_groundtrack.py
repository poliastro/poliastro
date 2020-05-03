import pytest

from poliastro.examples import iss
from poliastro.plotting.groundtrack import groundtrack


@pytest.mark.mpl_image_compare
def test_groundtrack_plotting():
    fig = groundtrack(iss, tof=5)

    return fig
