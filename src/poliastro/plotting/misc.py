import warnings

import astropy.units as u
import erfa
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
from poliastro.frames import Planes
from poliastro.plotting import OrbitPlotter


def plot_solar_system(
    epoch=None,
    labels=None,
    outer=True,
    backend=None,
    length_scale_units=u.km,
):
    """
    Plots the whole solar system in one single call.

    .. versionadded:: 0.9.0

    Parameters
    ----------
    epoch : ~astropy.time.Time, optional
        Epoch value of the plot, default to J2000.
    labels : list[str]
        A list of strings containing the labels of the bodies.
    outer : bool, optional
        Whether to print the outer Solar System, default to True.
    backend : ~poliastro.plotting.orbit.backends._base.OrbitPlotterBackend
        An instance of ``OrbitPlotterBackend`` for rendendering the scene.
    length_scale_units : ~astropy.units.Unit
        Desired units of lenght used for representing distances.

    Returns
    -------
    ~poliastro.plotting.orbit.plotter.OrbitPlotter
        An object for plotting orbits.

    """
    # Compute current epoch if none is provided
    if epoch is None:
        epoch = Time.now().tdb

    # Get a list of all bodies to be plotted in the scene
    bodies_list = [Mercury, Venus, Earth, Mars]
    if outer:
        bodies_list.extend([Jupiter, Saturn, Uranus, Neptune])

    # Assert same number of bodies and labels
    if labels is None:
        labels = [body.name for body in bodies_list]
    elif labels is not None and len(labels) != len(bodies_list):
        raise ValueError(f"A total of {len(bodies_list)} labels are required.")

    # Ignore warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", erfa.core.ErfaWarning)

        # Instantiate the plotter and set the desired reference frame
        plotter = OrbitPlotter(
            backend=backend,
            plane=Planes.EARTH_ECLIPTIC,
            length_scale_units=length_scale_units,
        )
        plotter.set_body_frame(Earth, epoch)

        # Plot desired Solar System bodies and return the plotter instance
        for body, label in zip(bodies_list, labels):
            plotter.plot_body_orbit(body, label=label, epoch=epoch)

    return plotter
