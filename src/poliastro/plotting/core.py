""" Plotting utilities.

"""
import warnings

import numpy as np
from astropy import units as u
from astropy.coordinates import CartesianRepresentation
from plotly.graph_objs import Layout, Scatter, Scatter3d, Surface

from poliastro.plotting.util import generate_circle, generate_sphere
from poliastro.util import norm

from ._base import BaseOrbitPlotter


class OrbitPlotter3D(BaseOrbitPlotter):
    """OrbitPlotter3D class.

    """

    def __init__(self, figure=None, dark=False):
        super().__init__(figure)
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

    def _plot_trajectory(self, trajectory, label, color, dashed):
        trace = Scatter3d(
            x=trajectory.x.to(u.km).value,
            y=trajectory.y.to(u.km).value,
            z=trajectory.z.to(u.km).value,
            name=label,
            line=dict(color=color, width=5, dash="dash" if dashed else "solid"),
            mode="lines",  # Boilerplate
        )
        self._figure.add_trace(trace)

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

        return self.figure


class OrbitPlotter2D(BaseOrbitPlotter):
    """OrbitPlotter2D class.

    .. versionadded:: 0.9.0
    """

    def __init__(self, figure=None):
        super().__init__(figure)
        self._layout = Layout(
            autosize=True,
            xaxis=dict(title="x (km)", constrain="domain"),
            yaxis=dict(title="y (km)", scaleanchor="x"),
            shapes=[],
        )

        self._frame = None

    def _project(self, rr):
        if self._frame is None:
            warnings.warn("Frame has not been set, using identity")
            frame = [[1, 0, 0] * u.one, [0, 1, 0] * u.one, [0, 0, 1] * u.one]
        else:
            frame = self._frame

        rr_proj = rr - rr.dot(frame[2])[:, None] * frame[2]
        x = rr_proj.dot(frame[0])
        y = rr_proj.dot(frame[1])
        return x, y

    def _plot_sphere(self, radius, color, name, center=[0, 0, 0] * u.km):
        xx, yy = generate_circle(radius, center)
        x_center, y_center = self._project(
            center[None]
        )  # Indexing trick to add one extra dimension

        # TODO: Review
        trace = Scatter(
            x=xx.to(u.km).value,
            y=yy.to(u.km).value,
            mode="markers",
            line=dict(color=color, width=5, dash="dash"),
            name=name,
            hoverinfo="none",  # TODO: Review
            showlegend=False,
        )

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
        self._figure.add_trace(trace)

    def _plot_trajectory(self, trajectory, label, color, dashed):
        rr = trajectory.represent_as(CartesianRepresentation).xyz.transpose()
        x, y = self._project(rr)

        trace = Scatter(
            x=x.to(u.km).value,
            y=y.to(u.km).value,
            name=label,
            line=dict(color=color, width=2, dash="dash" if dashed else "solid"),
            hoverinfo="none",  # TODO: Review
            mode="lines",  # Boilerplate
        )
        self._figure.add_trace(trace)

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
        if not self._frame:
            self.set_frame(*orbit.pqw())

        return super().plot(orbit, label=label, color=color)
