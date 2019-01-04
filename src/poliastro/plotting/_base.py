import os.path
from itertools import cycle

import numpy as np
import plotly.colors
from astropy import units as u
from astropy.coordinates import CartesianRepresentation
from plotly.offline import iplot, plot as export

from poliastro.plotting.util import BODY_COLORS, generate_label
from poliastro.util import norm


class BaseOrbitPlotter:
    """
    Parent Class for the 2D and 3D OrbitPlotter Classes based on Plotly.
    """

    def __init__(self):
        self._data = []  # type: List[dict]
        self._attractor = None
        self._attractor_data = {}  # type: dict
        self._attractor_radius = np.inf * u.km
        self._color_cycle = cycle(plotly.colors.DEFAULT_PLOTLY_COLORS)

    @property
    def figure(self):
        return dict(data=self._data + [self._attractor_data], layout=self._layout)

    def set_attractor(self, attractor):
        """Sets plotting attractor.

        Parameters
        ----------
        attractor : ~poliastro.bodies.Body
            Central body.
        """
        if self._attractor is None:
            self._attractor = attractor

        elif attractor is not self._attractor:
            raise NotImplementedError(
                "Attractor has already been set to {}.".format(self._attractor.name)
            )

    def _redraw_attractor(self, min_radius=0 * u.km):
        # Select a sensible value for the radius: realistic for low orbits,
        # visible for high and very high orbits
        radius = max(self._attractor.R.to(u.km), min_radius.to(u.km))
        if radius < self._attractor_radius:
            # If the resulting radius is smaller than the current one, redraw it
            shape = self._plot_sphere(
                radius,
                BODY_COLORS.get(self._attractor.name, "#999999"),
                self._attractor.name,
            )
            # Overwrite stored properties
            self._attractor_radius = radius
            self._attractor_data = shape

    def plot_trajectory(self, trajectory, *, label=None, color=None):
        """Plots a precomputed trajectory.

        An attractor must be set first.

        Parameters
        ----------
        trajectory : ~astropy.coordinates.CartesianRepresentation
            Trajectory to plot.
        label : string, optional
        color : string, optional

        """
        if self._attractor is None:
            raise ValueError(
                "An attractor must be set up first, please use "
                "set_attractor(Major_Body)."
            )
        else:
            self._redraw_attractor(
                trajectory.represent_as(CartesianRepresentation).norm().min() * 0.15
            )

        self._plot_trajectory(trajectory, str(label), color, False)

    def _plot_trajectory(self, trajectory, label, color, dashed):
        raise NotImplementedError

    def plot(self, orbit, *, label=None, color=None):
        """Plots state and osculating orbit in their plane.

        Parameters
        ----------
        orbit : ~poliastro.twobody.orbit.Orbit
        label : string, optional
        color : string, optional
        """
        if color is None:
            color = next(self._color_cycle)

        self.set_attractor(orbit.attractor)
        self._redraw_attractor(orbit.r_p * 0.15)  # Arbitrary threshold

        label = generate_label(orbit, label)
        trajectory = orbit.sample()

        self._plot_trajectory(trajectory, label, color, True)
        # Plot required 2D/3D shape in the position of the body
        radius = min(
            self._attractor_radius * 0.5, (norm(orbit.r) - orbit.attractor.R) * 0.3
        )  # Arbitrary thresholds
        shape = self._plot_sphere(radius, color, label, center=orbit.r)
        self._data.append(shape)

    def _prepare_plot(self, **layout_kwargs):
        # If there are no orbits, draw only the attractor
        if not self._data:
            self._redraw_attractor()

        if layout_kwargs:
            self._layout.update(layout_kwargs)

    def show(self, **layout_kwargs):
        """Shows the plot in the Notebook.
        """
        self._prepare_plot(**layout_kwargs)

        iplot(self.figure)

    def savefig(self, filename, **layout_kwargs):
        """Function for saving the plot locally.
        Parameters
        ----------
        filename : string, format = "anyname.anyformat"
                    anyformat can only be jpeg, png, svg, webp
        """
        self._prepare_plot(**layout_kwargs)

        basename, ext = os.path.splitext(filename)
        export(
            self.figure,
            image=ext[1:],
            image_filename=basename,
            show_link=False,  # Boilerplate
        )
