""" Plotting utilities.

"""
import os.path
from itertools import cycle

import numpy as np

from typing import List, Tuple  # flake8: noqa

import plotly.colors
from plotly.offline import iplot, plot as export
from plotly.graph_objs import Scatter3d, Surface, Layout, Scatter

from astropy import units as u
from astropy.coordinates import CartesianRepresentation

from poliastro.util import norm


BODY_COLORS = {
    "Sun": "#ffcc00",
    "Earth": "#204a87",
    "Jupiter": "#ba3821",
}


def _generate_label(orbit, label):
    orbit.epoch.out_subfmt = 'date_hm'
    label_ = '{}'.format(orbit.epoch.iso)
    if label:
        label_ += ' ({})'.format(label)

    return label_


def _generate_sphere(radius, center, num=20):
    u1 = np.linspace(0, 2 * np.pi, num)
    v1 = u1.copy()
    uu, vv = np.meshgrid(u1, v1)

    x_center, y_center, z_center = center

    xx = x_center + radius * np.cos(uu) * np.sin(vv)
    yy = y_center + radius * np.sin(uu) * np.sin(vv)
    zz = z_center + radius * np.cos(vv)

    return xx, yy, zz


class _BaseOrbitPlotter:
    """
    Parent Class for the 2D and 3D OrbitPlotter Classes based on Plotly.
    """

    def __init__(self):
        self._data = []   # type: List[dict]
        self._attractor = None
        self._attractor_data = {}  # type: dict
        self._attractor_radius = np.inf * u.km
        self._color_cycle = cycle(plotly.colors.DEFAULT_PLOTLY_COLORS)

    @property
    def figure(self):
        return dict(
            data=self._data + [self._attractor_data],
            layout=self._layout,
        )

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
            raise NotImplementedError("Attractor has already been set to {}.".format(self._attractor.name))

    def _redraw_attractor(self, min_radius=0 * u.km):
        # Select a sensible value for the radius: realistic for low orbits,
        # visible for high and very high orbits
        radius = max(self._attractor.R.to(u.km), min_radius.to(u.km))
        if radius < self._attractor_radius:
            # If the resulting radius is smaller than the current one, redraw it
            shape = self._plot_sphere(radius, BODY_COLORS.get(self._attractor.name, "#999999"), self._attractor.name)
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
            raise ValueError("An attractor must be set up first, please use "
                             "set_attractor(Major_Body).")
        else:
            self._redraw_attractor(trajectory.represent_as(CartesianRepresentation).norm().min() * 0.15)

        self._plot_trajectory(trajectory, str(label), color, False)

    def _plot_point(self, radius, color, name, center=[0, 0, 0] * u.km):
        return self._plot_sphere(radius, color, name, center)

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

        label = _generate_label(orbit, label)
        trajectory = orbit.sample()

        self._plot_trajectory(trajectory, label, color, True)
        # Plot required 2D/3D shape in the position of the body
        radius = min(self._attractor_radius * 0.5, (norm(orbit.r) - orbit.attractor.R) * 0.3)  # Arbitrary thresholds
        shape = self._plot_point(radius, color, label, center=orbit.r)
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
            image=ext[1:], image_filename=basename,
            show_link=False,  # Boilerplate
        )


class OrbitPlotter3D(_BaseOrbitPlotter):
    """OrbitPlotter3D class.
    """

    def __init__(self):
        super().__init__()
        self._layout = Layout(
            autosize=True,
            scene=dict(
                xaxis=dict(
                    title="x (km)",
                ),
                yaxis=dict(
                    title="y (km)",
                ),
                zaxis=dict(
                    title="z (km)",
                ),
                aspectmode="data",  # Important!
            ),
        )

    def _plot_sphere(self, radius, color, name, center=[0, 0, 0] * u.km):
        xx, yy, zz = _generate_sphere(radius, center)
        sphere = Surface(
            x=xx.to(u.km).value, y=yy.to(u.km).value, z=zz.to(u.km).value,
            name=name,
            colorscale=[[0, color], [1, color]],
            cauto=False, cmin=1, cmax=1, showscale=False,  # Boilerplate
        )
        return sphere

    @u.quantity_input(elev=u.rad, azim=u.rad, distance=u.km)
    def set_view(self, elev, azim, distance=5 * u.km):
        x = distance * np.cos(elev) * np.cos(azim)
        y = distance * np.cos(elev) * np.sin(azim)
        z = distance * np.sin(elev)

        self._layout.update({
            "scene": {
                "camera": {
                    "eye": {
                        "x": x.to(u.km).value, "y": y.to(u.km).value, "z": z.to(u.km).value
                    }
                }
            }
        })

    def _plot_trajectory(self, trajectory, label, color, dashed):
        trace = Scatter3d(
            x=trajectory.x.to(u.km).value, y=trajectory.y.to(u.km).value, z=trajectory.z.to(u.km).value,
            name=label,
            line=dict(
                color=color,
                width=5,
                dash='dash' if dashed else 'solid',
            ),
            mode="lines",  # Boilerplate
        )
        self._data.append(trace)


class OrbitPlotter2D(_BaseOrbitPlotter):
    """OrbitPlotter2D class.

    .. versionadded:: 0.9.0
    """

    def __init__(self):
        super().__init__()
        self._layout = Layout(
            autosize=True,
            xaxis=dict(
                title="x (km)",
                constrain="domain",
            ),
            yaxis=dict(
                title="y (km)",
                scaleanchor="x",
            ),
        )
        self._layout.update({
            "shapes": []
        })
        self._frame = None

    def _plot_point(self, radius, color, name, center=[0, 0, 0] * u.km):
        x_center, y_center = self._project(center[None])  # Indexing trick to add one extra dimension

        trace = Scatter(x=x_center.to(u.km).value, y=y_center.to(u.km).value, mode='markers',
                        marker=dict(size=10, color=color), name=name)
        return trace

    def _plot_sphere(self, radius, color, name, center=[0, 0, 0] * u.km):
        x_center, y_center = self._project(center[None])  # Indexing trick to add one extra dimension

        trace = {
            'type': 'circle',
            'xref': 'x',
            'yref': 'y',
            'x0': (x_center[0] - radius).to(u.km).value,
            'y0': (y_center[0] - radius).to(u.km).value,
            'x1': (x_center[0] + radius).to(u.km).value,
            'y1': (y_center[0] + radius).to(u.km).value,
            'opacity': 1,
            'fillcolor': color,
            'line': {
                'color': color,
            },
        }

        self._layout["shapes"] += (trace,)
        return {}  # This stores {} in self._data already

    def _plot_trajectory(self, trajectory, label, color, dashed):
        rr = trajectory.represent_as(CartesianRepresentation).xyz.transpose()
        x, y = self._project(rr)

        trace = Scatter(
            x=x.to(u.km).value, y=y.to(u.km).value,
            name=label,
            line=dict(
                color=color,
                width=2,
                dash='dash' if dashed else 'solid',
            ),
            mode="lines",  # Boilerplate
        )
        self._data.append(trace)

    def _project(self, rr):
        rr_proj = rr - rr.dot(self._frame[2])[:, None] * self._frame[2]
        x = rr_proj.dot(self._frame[0])
        y = rr_proj.dot(self._frame[1])
        return x, y

    def set_frame(self, p_vec, q_vec, w_vec):
        """Sets perifocal frame.

        Raises
        ------
        ValueError
            If the vectors are not a set of mutually orthogonal unit vectors.
        """
        if self._frame and self._data:
            raise NotImplementedError("OrbitPlotter2D does not support reprojecting yet")

        if not np.allclose([norm(v) for v in (p_vec, q_vec, w_vec)], 1):
            raise ValueError("Vectors must be unit.")
        elif not np.allclose([p_vec.dot(q_vec),
                              q_vec.dot(w_vec),
                              w_vec.dot(p_vec)], 0):
            raise ValueError("Vectors must be mutually orthogonal.")
        else:
            self._frame = p_vec, q_vec, w_vec

    def plot(self, orbit, *, label=None, color=None):
        if not self._frame:
            self.set_frame(*orbit.pqw())

        super().plot(orbit, label=label, color=color)


def plot(state, label=None, color=None, plotter=OrbitPlotter2D):
    """Quickly plots an :py:class:`~poliastro.twobody.orbit.Orbit`.

    For more advanced tuning, use the :py:class:`OrbitPlotter2D` class
    and similar ones.

    """
    op = plotter()
    op.plot(state, label=label, color=color)
    op.show()

    return op


def plot3d(orbit, *, label=None, color=None):
    """Plots an :py:class:`~poliastro.twobody.orbit.Orbit` in 3D.

    For more advanced tuning, use the :py:class:`OrbitPlotter3D` class.

    """
    return plot(orbit, label=label, color=color, plotter=OrbitPlotter3D)
