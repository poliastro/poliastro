import matplotlib

from poliastro.examples import iss
from poliastro.plotting.static import StaticOrbitPlotter


def test_type_of_anim():
    op = StaticOrbitPlotter()
    ss = iss
    k = op.anim(ss)
    assert type(k) == matplotlib.animation.FuncAnimation
