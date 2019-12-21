from typing import Union

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

from .core import OrbitPlotter2D, OrbitPlotter3D
from .static import StaticOrbitPlotter


def _plot_bodies(orbit_plotter, outer=True, epoch=None):
    bodies = [Mercury, Venus, Earth, Mars]
    if outer:
        bodies.extend([Jupiter, Saturn, Uranus, Neptune])

    for body in bodies:
        orbit_plotter.plot(body.get_mean_orbit(epoch), label=str(body))


def _plot_solar_system_2d(outer=True, epoch=None, interactive=False):
    pqw = Earth.get_mean_orbit().pqw()
    if interactive:
        orbit_plotter = (
            OrbitPlotter2D()
        )  # type: Union[OrbitPlotter2D, StaticOrbitPlotter]
        orbit_plotter.set_frame(*pqw)
    else:
        orbit_plotter = StaticOrbitPlotter()
        orbit_plotter.set_frame(*pqw)

    _plot_bodies(orbit_plotter, outer, epoch)

    return orbit_plotter


def _plot_solar_system_3d(outer=True, epoch=None):
    orbit_plotter = OrbitPlotter3D()

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
    if not interactive and use_3d:
        raise ValueError(
            "The static plotter does not support 3D, use `interactive=True`"
        )
    elif use_3d:
        op = _plot_solar_system_3d(outer, epoch)
    else:
        op = _plot_solar_system_2d(outer, epoch, interactive)

    return op
