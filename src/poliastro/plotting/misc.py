from typing import Union

from astropy.time import Time

from poliastro.bodies import (
    Earth,
    Jupiter,
    Mars,
    Mercury,
    Neptune,
    Saturn,
    Uranus,
    Venus,
)

from ..frames import Planes
from .core import OrbitPlotter2D, OrbitPlotter3D
from .static import StaticOrbitPlotter


def _plot_bodies(orbit_plotter, outer=True, epoch=None):
    bodies = [Mercury, Venus, Earth, Mars]
    if outer:
        bodies.extend([Jupiter, Saturn, Uranus, Neptune])

    for body in bodies:
        orbit_plotter.plot_body_orbit(body, epoch)


def _plot_solar_system_2d(epoch, outer=True, interactive=False):
    if interactive:
        orbit_plotter = OrbitPlotter2D(
            plane=Planes.EARTH_ECLIPTIC
        )  # type: Union[OrbitPlotter2D, StaticOrbitPlotter]
    else:
        orbit_plotter = StaticOrbitPlotter(plane=Planes.EARTH_ECLIPTIC)

    orbit_plotter.set_body_frame(Earth, epoch)

    _plot_bodies(orbit_plotter, outer, epoch)

    return orbit_plotter


def _plot_solar_system_3d(epoch, outer=True):
    orbit_plotter = OrbitPlotter3D(plane=Planes.EARTH_ECLIPTIC)

    _plot_bodies(orbit_plotter, outer, epoch)

    return orbit_plotter


def plot_solar_system(outer=True, epoch=None, use_3d=False, interactive=False):
    """
    Plots the whole solar system in one single call.

    .. versionadded:: 0.9.0

    Parameters
    ------------
    outer : bool, optional
        Whether to print the outer Solar System, default to True.
    epoch : ~astropy.time.Time, optional
        Epoch value of the plot, default to J2000.
    use_3d : bool, optional
        Produce 3D plot, default to False.
    interactive : bool, optional
        Produce an interactive (rather than static) image of the orbit, default to False.
        This option requires Plotly properly installed and configured for your environment.

    """
    if epoch is None:
        epoch = Time.now().tdb

    if not interactive and use_3d:
        raise ValueError(
            "The static plotter does not support 3D, use `interactive=True`"
        )
    elif use_3d:
        op = _plot_solar_system_3d(epoch, outer)
    else:
        op = _plot_solar_system_2d(epoch, outer, interactive)

    return op
