""" Plotting utilities.

"""
from itertools import cycle

import numpy as np
import plotly.colors
from astropy import units as u
from astropy.coordinates import CartesianRepresentation
from plotly.graph_objects import Figure, Layout, Scatter, Scatter3d, Surface

from poliastro.plotting.util import generate_sphere
from poliastro.util import norm

from ._base import BaseOrbitPlotter


class _PlotlyOrbitPlotter(BaseOrbitPlotter):
    def __init__(self, figure=None, *, num_points=150):
        super().__init__(num_points)

        self._figure = figure or Figure()
        self._layout = None

        self._color_cycle = cycle(plotly.colors.DEFAULT_PLOTLY_COLORS)

    def _prepare_plot(self):
        super()._prepare_plot()

        self._figure.layout.update(self._layout)

    def _get_colors(self, color):
        if color is None:
            color = next(self._color_cycle)

        return [color]

    def plot_trajectory(self, trajectory, *, label=None, color=None):
        super().plot_trajectory(trajectory, label=label, color=color)

        if not self._figure._in_batch_mode:
            return self.show()

    def plot(self, orbit, *, label=None, color=None):
        super().plot(orbit, label=label, color=color)

        if not self._figure._in_batch_mode:
            return self.show()

    def show(self):
        """Shows the plot in the Notebook.

        Updates the layout and returns the underlying figure.

        """
        self._prepare_plot()
        return self._figure


class OrbitPlotter3D(_PlotlyOrbitPlotter):
    """OrbitPlotter3D class.

    """

    def __init__(self, figure=None, dark=False, *, num_points=150):
        super().__init__(figure, num_points=num_points)
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

    def _plot_point(self, radius, color, name, center=[0, 0, 0] * u.km):
        # We use _plot_sphere here because it's not easy to specify the size of a marker
        # in data units instead of pixels, see
        # https://stackoverflow.com/q/47086547
        return self._plot_sphere(radius, color, name, center)

    def _plot_sphere(self, radius, color, name, center=[0, 0, 0] * u.km):
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

    def _plot_trajectory(self, trajectory, label, colors, dashed):
        trace = Scatter3d(
            x=trajectory.x.to(u.km).value,
            y=trajectory.y.to(u.km).value,
            z=trajectory.z.to(u.km).value,
            name=label,
            line=dict(color=colors[0], width=5, dash="dash" if dashed else "solid"),
            mode="lines",  # Boilerplate
        )
        self._figure.add_trace(trace)

        return trace, [trace.line.color]

    @u.quantity_input(elev=u.rad, azim=u.rad, distance=u.km)
    def set_view(self, elev, azim, distance=5 * u.km):
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


class OrbitPlotter2D(_PlotlyOrbitPlotter):
    """OrbitPlotter2D class.

    .. versionadded:: 0.9.0
    """

    def __init__(self, figure=None, *, num_points=150):
        super().__init__(figure, num_points=num_points)
        self._layout = Layout(
            autosize=True,
            xaxis=dict(title="x (km)", constrain="domain"),
            yaxis=dict(title="y (km)", scaleanchor="x"),
            shapes=[],
        )

        self._frame = None

    def _project(self, rr):
        rr_proj = rr - rr.dot(self._frame[2])[:, None] * self._frame[2]
        x = rr_proj.dot(self._frame[0])
        y = rr_proj.dot(self._frame[1])
        return x, y

    def _plot_point(self, radius, color, name, center=[0, 0, 0] * u.km):
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

    def _plot_sphere(self, radius, color, name, center=[0, 0, 0] * u.km):
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

    def _plot_trajectory(self, trajectory, label, colors, dashed):
        if self._frame is None:
            raise ValueError(
                "A frame must be set up first, please use "
                "set_frame(*orbit.pqw()) or plot(orbit)"
            )

        rr = trajectory.represent_as(CartesianRepresentation).xyz.transpose()
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

        return trace, [trace.line.color]

    def set_frame(self, p_vec, q_vec, w_vec):
        """Sets perifocal frame.

        Raises
        ------
        ValueError
            If the vectors are not a set of mutually orthogonal unit vectors.

        """
        if self._frame and self.trajectories:
            raise NotImplementedError(
                "OrbitPlotter2D does not support reprojecting yet"
            )

        if not np.allclose([norm(v) for v in (p_vec, q_vec, w_vec)], 1):
            raise ValueError("Vectors must be unit.")
        elif not np.allclose([p_vec.dot(q_vec), q_vec.dot(w_vec), w_vec.dot(p_vec)], 0):
            raise ValueError("Vectors must be mutually orthogonal.")
        else:
            self._frame = p_vec, q_vec, w_vec

    def plot(self, orbit, *, label=None, color=None):
        """Plots state and osculating orbit in their plane.

        Parameters
        ----------
        orbit : ~poliastro.twobody.orbit.Orbit
            Orbit to plot.
        label : string, optional
            Label of the orbit.
        color : string, optional
            Color of the line and the position.

        """
        if not self._frame:
            self.set_frame(*orbit.pqw())

        return super().plot(orbit, label=label, color=color)
