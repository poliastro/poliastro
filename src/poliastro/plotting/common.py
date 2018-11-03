""" Plotting utilities.

"""
from typing import List, Tuple  # flake8: noqa

from poliastro.bodies import (Earth, Jupiter, Mars, Mercury, Neptune, Saturn,
                              Uranus, Venus)
from poliastro.twobody.orbit import Orbit

from .interactive import OrbitPlotter2D, OrbitPlotter3D


BODY_COLORS = {
    "Sun": "#ffcc00",
    "Earth": "#204a87",
    "Jupiter": "#ba3821",
}


def plot_solar_system(outer=True, epoch=None, plotter=OrbitPlotter2D):
    """
    Plots the whole solar system in one single call.

    .. versionadded:: 0.9.0

    Parameters
    ------------
    outer : bool, optional
        Whether to print the outer Solar System, default to True.
    epoch: ~astropy.time.Time, optional
        Epoch value of the plot, default to J2000.
    """
    bodies = [Mercury, Venus, Earth, Mars]
    if outer:
        bodies.extend([Jupiter, Saturn, Uranus, Neptune])

    op = plotter()
    for body in bodies:
        orb = Orbit.from_body_ephem(body, epoch)
        op.plot(orb, label=str(body))

    # Sets frame to the orbit of the Earth by default
    # FIXME: https://github.com/poliastro/poliastro/issues/316
    op.set_frame(*Orbit.from_body_ephem(Earth, epoch).pqw())

    return op
