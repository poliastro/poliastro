import matplotlib

from poliastro.examples import iss
from poliastro.plotting import AnimatedOrbitPlotter


def test_type_of_anim():
    op = AnimatedOrbitPlotter()
    ss = iss
    k = op.anim(ss)
    assert type(k) == matplotlib.animation.FuncAnimation
