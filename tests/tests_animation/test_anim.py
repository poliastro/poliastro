import matplotlib
from astropy.time import Time
from matplotlib import pyplot as plt

from poliastro.examples import churi, iss, molniya
from poliastro.frames import Planes
from poliastro.plotting.static import StaticOrbitPlotter


def test_type_of_anim():
    op = StaticOrbitPlotter()
    ss = iss
    k = op.anim(ss)
    assert type(k) == matplotlib.animation.FuncAnimation
