""" Plotting utilities.

"""
from itertools import cycle

import numpy as np
import plotly.colors
from astropy import units as u
from plotly.graph_objects import Figure, Layout, Scatter, Scatter3d, Surface

from poliastro.plotting.util import generate_sphere

from ..frames import Planes
from ._base import BaseOrbitPlotter, Mixin2D


class _PlotlyOrbitPlotter(BaseOrbitPlotter):
    def __init__(self, figure=None, *, num_points=150, plane=None):
        super().__init__(num_points=num_points, plane=plane)

        self._figure = figure or Figure()
        self._layout = None

        self._color_cycle = cycle(plotly.colors.DEFAULT_PLOTLY_COLORS)

    def _clear_attractor(self):
        # FIXME: Implement
        pass

    def _get_colors(self, color, trail):
        # TODO: Support trail
        if color is None:
            color = next(self._color_cycle)

        return [color]

    def plot_trajectory(self, coordinates, *, label=None, color=None, trail=False):
        """Plots a precomputed trajectory.

        An attractor must be set first.

        Parameters
        ----------
        coordinates : ~astropy.coordinates.CartesianRepresentation
            Trajectory to plot.
        label : str, optional
            Label of the trajectory.
        color : str, optional
            Color of the trajectory.
        trail : bool, optional
            Fade the orbit trail, default to False.

        """

        super().plot_trajectory(coordinates, label=label, color=color, trail=trail)

        if not self._figure._in_batch_mode:
            return self.show()

    def plot(self, orbit, *, label=None, color=None, trail=False):
        """Plots state and osculating orbit in their plane.

        Parameters
        ----------
        orbit : ~poliastro.twobody.orbit.Orbit
            Orbit to plot.
        label : str, optional
            Label of the orbit.
        color : str, optional
            Color of the line and the position.
        trail : bool, optional
            Fade the orbit trail, default to False.

        """
        super().plot(orbit, label=label, color=color, trail=trail)

        if not self._figure._in_batch_mode:
            return self.show()

    def plot_body_orbit(
        self,
        body,
        epoch,
        *,
        label=None,
        color=None,
        trail=False,
    ):
        """Plots complete revolution of body and current position.

        Parameters
        ----------
        body : poliastro.bodies.SolarSystemPlanet
            Body.
        epoch : astropy.time.Time
            Epoch of current position.
        label : str, optional
            Label of the orbit, default to the name of the body.
        color : str, optional
            Color of the line and the position.
        trail : bool, optional
            Fade the orbit trail, default to False.

        """
        super().plot_body_orbit(body, epoch, label=label, color=color, trail=trail)

        if not self._figure._in_batch_mode:
            return self.show()

    def plot_ephem(self, ephem, epoch=None, *, label=None, color=None, trail=False):
        """Plots Ephem object over its sampling period.

        Parameters
        ----------
        ephem : ~poliastro.ephem.Ephem
            Ephemerides to plot.
        epoch : astropy.time.Time, optional
            Epoch of the current position, none will be used if not given.
        label : str, optional
            Label of the orbit, default to the name of the body.
        color : str, optional
            Color of the line and the position.
        trail : bool, optional
            Fade the orbit trail, default to False.

        """

        super().plot_ephem(ephem, epoch, label=label, color=color, trail=trail)

        if not self._figure._in_batch_mode:
            return self.show()

    def show(self):
        """Shows the plot in the Notebook.

        Updates the layout and returns the underlying figure.

        """
        if self._attractor is not None:
            self._redraw_attractor()

        self._figure.layout.update(self._layout)
        return self._figure


class OrbitPlotter3D(_PlotlyOrbitPlotter):
    """OrbitPlotter3D class."""

    def __init__(self, figure=None, dark=False, *, num_points=150, plane=None):
        super().__init__(figure, num_points=num_points, plane=plane)
        self._layout = Layout(
            autosize=True,
            scene=dict(
                xaxis=dict(title="x (km)"),
                yaxis=dict(title="y (km)"),
                zaxis=dict(title="z (km)"),
                aspectmode="data",  # Important!
            ),
        )
        if dark:
            self._layout.template = "plotly_dark"

    def _draw_point(self, radius, color, name, center=[0, 0, 0] * u.km):
        # We use _plot_sphere here because it's not easy to specify the size of a marker
        # in data units instead of pixels, see
        # https://stackoverflow.com/q/47086547
        return self._draw_sphere(radius, color, name, center)

    def _draw_sphere(self, radius, color, name, center=[0, 0, 0] * u.km):
        xx, yy, zz = generate_sphere(radius, center)
        sphere = Surface(
            x=xx.to(u.km).value,
            y=yy.to(u.km).value,
            z=zz.to(u.km).value,
            name=name,
            colorscale=[[0, color], [1, color]],
            cauto=False,
            cmin=1,
            cmax=1,
            showscale=False,
        )
        self._figure.add_trace(sphere)

        return sphere

    def _plot_coordinates(self, coordinates, label, colors, dashed):
        trace = Scatter3d(
            x=coordinates.x.to(u.km).value,
            y=coordinates.y.to(u.km).value,
            z=coordinates.z.to(u.km).value,
            name=label,
            line=dict(color=colors[0], width=5, dash="dash" if dashed else "solid"),
            mode="lines",  # Boilerplate
        )
        self._figure.add_trace(trace)

        return trace, [trace.line.color]

    def plot(self, orbit, *, label=None, color=None, trail=False):
        """Plots state and osculating orbit in their plane.

        Parameters
        ----------
        orbit : ~poliastro.twobody.orbit.Orbit
            Orbit to plot.
        label : str, optional
            Label of the orbit.
        color : str, optional
            Color of the line and the position.
        trail : bool, optional
            Fade the orbit trail, default to False.

        """
        if trail:
            raise NotImplementedError("trail not supported yet")

        return super().plot(orbit, label=label, color=color, trail=trail)

    @u.quantity_input(elev=u.rad, azim=u.rad, distance=u.km)
    def set_view(self, elev, azim, distance=5 * u.km):
        """Changes 3D view."""
        x = distance * np.cos(elev) * np.cos(azim)
        y = distance * np.cos(elev) * np.sin(azim)
        z = distance * np.sin(elev)

        self._layout.update(
            {
                "scene": {
                    "camera": {
                        "eye": {
                            "x": x.to(u.km).value,
                            "y": y.to(u.km).value,
                            "z": z.to(u.km).value,
                        }
                    }
                }
            }
        )

        if not self._figure._in_batch_mode:
            return self.show()


class OrbitPlotter2D(_PlotlyOrbitPlotter, Mixin2D):
    """OrbitPlotter2D class.

    .. versionadded:: 0.9.0
    """

    def __init__(self, figure=None, *, num_points=150, plane=None):
        super().__init__(figure, num_points=num_points, plane=plane)
        self._layout = Layout(
            autosize=True,
            xaxis=dict(title="x (km)", constrain="domain"),
            yaxis=dict(title="y (km)", scaleanchor="x"),
            shapes=[],
        )

        self._frame = None

    def _redraw(self):
        raise NotImplementedError("OrbitPlotter2D does not support reprojecting yet")

    def _draw_point(self, radius, color, name, center=[0, 0, 0] * u.km):
        x_center, y_center = self._project(
            center[None]
        )  # Indexing trick to add one extra dimension

        trace = Scatter(
            x=x_center.to(u.km).value,
            y=y_center.to(u.km).value,
            mode="markers",
            marker=dict(size=10, color=color),
            name=name,
        )
        self._figure.add_trace(trace)

        return trace

    def _draw_sphere(self, radius, color, name, center=[0, 0, 0] * u.km):
        x_center, y_center = self._project(
            center[None]
        )  # Indexing trick to add one extra dimension

        shape = {
            "type": "circle",
            "xref": "x",
            "yref": "y",
            "x0": (x_center[0] - radius).to(u.km).value,
            "y0": (y_center[0] - radius).to(u.km).value,
            "x1": (x_center[0] + radius).to(u.km).value,
            "y1": (y_center[0] + radius).to(u.km).value,
            "opacity": 1,
            "fillcolor": color,
            "line": {"color": color},
        }

        self._layout.shapes += (shape,)

        return shape

    def _plot_coordinates(self, coordinates, label, colors, dashed):
        if self._frame is None:
            raise ValueError(
                "A frame must be set up first, please use "
                "set_orbit_frame(orbit) or plot(orbit)"
            )

        rr = coordinates.xyz.transpose()
        x, y = self._project(rr)

        trace = Scatter(
            x=x.to(u.km).value,
            y=y.to(u.km).value,
            name=label,
            line=dict(color=colors[0], width=2, dash="dash" if dashed else "solid"),
            hoverinfo="none",  # TODO: Review
            mode="lines",  # Boilerplate
        )
        self._figure.add_trace(trace)

        return trace

    def plot_trajectory(self, coordinates, *, label=None, color=None, trail=False):
        """Plots a precomputed trajectory.

        An attractor must be set first.

        Parameters
        ----------
        coordinates : ~astropy.coordinates.CartesianRepresentation
            Trajectory to plot.
        label : str, optional
            Label of the trajectory.
        color : str, optional
            Color of the trajectory.
        trail : bool, optional
            Fade the orbit trail, default to False.

        """
        if self._frame is None:
            raise ValueError(
                "A frame must be set up first, please use "
                "set_orbit_frame(orbit) or plot(orbit)"
            )

        return super().plot_trajectory(
            coordinates, label=label, color=color, trail=trail
        )

    def plot_ephem(self, ephem, epoch=None, *, label=None, color=None, trail=False):
        """Plots Ephem object over its sampling period.

        Parameters
        ----------
        ephem : ~poliastro.ephem.Ephem
            Ephemerides to plot.
        epoch : astropy.time.Time, optional
            Epoch of the current position, none will be used if not given.
        label : str, optional
            Label of the orbit, default to the name of the body.
        color : str, optional
            Color of the line and the position.
        trail : bool, optional
            Fade the orbit trail, default to False.

        """
        if self._frame is None:
            raise ValueError(
                "A frame must be set up first, please use "
                "set_orbit_frame(orbit) or plot(orbit)"
            )
        super().plot_ephem(ephem, epoch, label=label, color=color, trail=trail)

    def plot(self, orbit, *, label=None, color=None, trail=False):
        """Plots state and osculating orbit in their plane.

        Parameters
        ----------
        orbit : ~poliastro.twobody.orbit.Orbit
            Orbit to plot.
        label : str, optional
            Label of the orbit.
        color : str, optional
            Color of the line and the position.
        trail : bool, optional
            Fade the orbit trail, default to False.

        """

        if trail:
            raise NotImplementedError("trail not supported yet")

        if not self._frame:
            self.set_orbit_frame(orbit)

        return super().plot(orbit, label=label, color=color, trail=trail)

    def plot_body_orbit(
        self,
        body,
        epoch,
        plane=Planes.EARTH_ECLIPTIC,
        *,
        label=None,
        color=None,
        trail=False,
    ):
        """Plots complete revolution of body and current position.

        Parameters
        ----------
        body : poliastro.bodies.SolarSystemPlanet
            Body.
        epoch : astropy.time.Time
            Epoch of current position.
        plane : ~poliastro.frames.enums.Planes
            Reference plane.
        label : str, optional
            Label of the orbit, default to the name of the body.
        color : str, optional
            Color of the line and the position.
        trail : bool, optional
            Fade the orbit trail, default to False.

        """
        if self._frame is None:
            self.set_body_frame(body, epoch)

        return super().plot_body_orbit(
            body, epoch, label=label, color=color, trail=trail
        )
