""" Plotting utilities.

"""
import numpy as np
from astropy import units as u
from plotly.graph_objs import Layout, Scatter, Scatter3d, Surface

from poliastro.plotting.util import generate_circle, generate_sphere

from ._base import BaseOrbitPlotter


class OrbitPlotter3D(BaseOrbitPlotter):
    """OrbitPlotter3D class.
    """

    def __init__(self, dark=False):
        super().__init__()
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
            showscale=False,  # Boilerplate
        )
        return sphere

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

    def _plot_trajectory(self, trajectory, label, color, dashed):
        trace = Scatter3d(
            x=trajectory.x.to(u.km).value,
            y=trajectory.y.to(u.km).value,
            z=trajectory.z.to(u.km).value,
            name=label,
            line=dict(color=color, width=5, dash="dash" if dashed else "solid"),
            mode="lines",  # Boilerplate
        )
        self._data.append(trace)


class OrbitPlotter2D(BaseOrbitPlotter):
    """OrbitPlotter2D class.

    .. versionadded:: 0.9.0
    """

    def __init__(self):
        super().__init__()
        self._layout = Layout(
            autosize=True,
            xaxis=dict(title="x (km)", constrain="domain"),
            yaxis=dict(title="y (km)", scaleanchor="x"),
        )
        self._layout.update({"shapes": []})

    def _plot_sphere(self, radius, color, name, center=[0, 0, 0] * u.km):
        xx, yy = generate_circle(radius, center)
        x_center, y_center, z_center = center
        trace = Scatter(
            x=xx.to(u.km).value,
            y=yy.to(u.km).value,
            mode="markers",
            line=dict(color=color, width=5, dash="dash"),
            name=name,
        )
        self._layout["shapes"] += (
            {
                "type": "circle",
                "xref": "x",
                "yref": "y",
                "x0": (x_center - radius).to(u.km).value,
                "y0": (y_center - radius).to(u.km).value,
                "x1": (x_center + radius).to(u.km).value,
                "y1": (y_center + radius).to(u.km).value,
                "opacity": 1,
                "fillcolor": color,
                "line": {"color": color},
            },
        )
        return trace

    def _plot_trajectory(self, trajectory, label, color, dashed):
        trace = Scatter(
            x=trajectory.x.to(u.km).value,
            y=trajectory.y.to(u.km).value,
            name=label,
            line=dict(color=color, width=2, dash="dash" if dashed else "solid"),
            mode="lines",  # Boilerplate
        )
        self._data.append(trace)
