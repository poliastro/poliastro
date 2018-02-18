""" Plotting utilities.

"""
import os.path
from itertools import cycle

import numpy as np

from typing import List, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt

import plotly.colors
from plotly.offline import iplot, plot as export
from plotly.graph_objs import Scatter3d, Surface, Layout

from astropy import units as u

from poliastro.util import norm
from poliastro.bodies import (Earth, Jupiter, Mars, Mercury, Neptune, Saturn,
                              Uranus, Venus)
from poliastro.twobody.orbit import Orbit


BODY_COLORS = {
    "Sun": "#ffcc00",
    "Earth": "#204a87",
    "Jupiter": "#ba3821",
}


def plot(state, label=None, color=None):
    """Plots a ``State``.

    For more advanced tuning, use the :py:class:`OrbitPlotter` class.

    """
    op = OrbitPlotter()
    return op.plot(state, label=label, color=color)


def plot3d(orbit, *, label=None, color=None):
    frame = OrbitPlotter3D()
    frame.plot(orbit, label=label, color=color)
    frame.show()

    return frame


class OrbitPlotter(object):
    """OrbitPlotter class.

    This class holds the perifocal plane of the first
    :py:class:`~poliastro.twobody.State` plotted in it using
    :py:meth:`plot`, so all following
    plots will be projected on that plane. Alternatively, you can call
    :py:meth:`set_frame` to set the frame before plotting.

    """

    def __init__(self, ax=None, num_points=100):
        """Constructor.

        Parameters
        ----------
        ax : ~matplotlib.axes.Axes
            Axes in which to plot. If not given, new ones will be created.
        num_points : int, optional
            Number of points to use in plots, default to 100.

        """
        self.ax = ax
        if not self.ax:
            _, self.ax = plt.subplots(figsize=(6, 6))
        self.num_points = num_points
        self._frame = None
        self._attractor_radius = None
        self._orbits = list(tuple())  # type: List[Tuple[Orbit, str]]

    @property
    def orbits(self):
        return self._orbits

    def set_frame(self, p_vec, q_vec, w_vec):
        """Sets perifocal frame.

        Raises
        ------
        ValueError
            If the vectors are not a set of mutually orthogonal unit vectors.

        """
        if not np.allclose([norm(v) for v in (p_vec, q_vec, w_vec)], 1):
            raise ValueError("Vectors must be unit.")
        elif not np.allclose([p_vec.dot(q_vec),
                              q_vec.dot(w_vec),
                              w_vec.dot(p_vec)], 0):
            raise ValueError("Vectors must be mutually orthogonal.")
        else:
            self._frame = p_vec, q_vec, w_vec

        if self._orbits:
            self._redraw()

    def _redraw(self):
        for artist in self.ax.lines + self.ax.collections:
            artist.remove()
        self._attractor_radius = None
        for orbit, label in self._orbits:
            self.plot(orbit, label)
        self.ax.relim()
        self.ax.autoscale()

    def set_attractor(self, orbit):
        """Sets plotting attractor.

        Parameters
        ----------
        orbit : ~poliastro.twobody.orbit.Orbit
            orbit with attractor to plot.

        """
        radius = max(orbit.attractor.R.to(u.km).value,
                     orbit.r_p.to(u.km).value / 6)
        color = BODY_COLORS.get(orbit.attractor.name, "#999999")
        self._attractor_radius = radius

        for attractor in self.ax.findobj(match=mpl.patches.Circle):
            attractor.remove()

        self.ax.add_patch(mpl.patches.Circle((0, 0), radius, lw=0, color=color))

    def plot(self, orbit, label=None, color=None):
        """Plots state and osculating orbit in their plane.

        """
        if not self._frame:
            self.set_frame(*orbit.pqw())

        if (orbit, label) not in self._orbits:
            self._orbits.append((orbit, label))

        # if new attractor radius is smaller, plot it
        new_radius = max(orbit.attractor.R.to(u.km).value,
                         orbit.r_p.to(u.km).value / 6)
        if not self._attractor_radius:
            self.set_attractor(orbit)
        elif new_radius < self._attractor_radius:
            self.set_attractor(orbit)

        lines = []

        _, positions = orbit.sample(self.num_points)
        rr = positions.get_xyz().transpose()

        # Project on OrbitPlotter frame
        # x_vec, y_vec, z_vec = self._frame
        rr_proj = rr - rr.dot(self._frame[2])[:, None] * self._frame[2]
        x = rr_proj.dot(self._frame[0])
        y = rr_proj.dot(self._frame[1])

        # Plot current position
        l, = self.ax.plot(x[0].to(u.km).value, y[0].to(u.km).value,
                          'o', mew=0, color=color)
        lines.append(l)

        l, = self.ax.plot(x.to(u.km).value, y.to(u.km).value,
                          '--', color=l.get_color())
        lines.append(l)

        if label:
            # This will apply the label to either the point or the osculating
            # orbit depending on the last plotted line, as they share variable
            if not self.ax.get_legend():
                size = self.ax.figure.get_size_inches() + [8, 0]
                self.ax.figure.set_size_inches(size)
            orbit.epoch.out_subfmt = 'date_hm'
            label = '{} ({})'.format(orbit.epoch.iso, label)

            l.set_label(label)
            self.ax.legend(bbox_to_anchor=(1.05, 1), title="Names and epochs")

        self.ax.set_xlabel("$x$ (km)")
        self.ax.set_ylabel("$y$ (km)")
        self.ax.set_aspect(1)

        return lines


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


def _plot_sphere(radius, color, name, center=[0, 0, 0] * u.km):
    xx, yy, zz = _generate_sphere(radius, center)
    sphere = Surface(
        x=xx.to(u.km).value, y=yy.to(u.km).value, z=zz.to(u.km).value,
        name=name,
        colorscale=[[0, color], [1, color]],
        cauto=False, cmin=1, cmax=1, showscale=False,  # Boilerplate
    )
    return sphere


class OrbitPlotter3D:
    """OrbitPlotter3D class.

    """

    def __init__(self):
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
        self._data = []  # type: List[dict]

        # TODO: Refactor?
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

    def _redraw_attractor(self, min_radius=0 * u.km):
        # Select a sensible value for the radius: realistic for low orbits,
        # visible for high and very high orbits
        radius = max(self._attractor.R.to(u.km), min_radius.to(u.km))

        # If the resulting radius is smaller than the current one, redraw it
        if radius < self._attractor_radius:
            sphere = _plot_sphere(
                radius, BODY_COLORS.get(self._attractor.name, "#999999"), self._attractor.name)

            # Overwrite stored properties
            self._attractor_radius = radius
            self._attractor_data = sphere

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

    @u.quantity_input(elev=u.rad, azim=u.rad)
    def set_view(self, elev, azim, distance=5):
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

    def plot(self, orbit, *, label=None, color=None):
        """Plots state and osculating orbit in their plane.

        """
        if color is None:
            color = next(self._color_cycle)

        self.set_attractor(orbit.attractor)

        self._redraw_attractor(orbit.r_p * 0.15)  # Arbitrary threshold

        label = _generate_label(orbit, label)
        _, trajectory = orbit.sample()

        self._plot_trajectory(trajectory, label, color, True)

        # Plot sphere in the position of the body
        radius = min(self._attractor_radius * 0.5, (norm(orbit.r) - orbit.attractor.R) * 0.3)  # Arbitrary thresholds
        sphere = _plot_sphere(
            radius, color, label, center=orbit.r)

        self._data.append(sphere)

    def plot_trajectory(self, trajectory, *, label=None, color=None):
        """Plots a precomputed trajectory.

        An attractor must be set first.

        Parameters
        ----------
        trajectory : ~astropy.coordinates.CartesianRepresentation
            Trajectory to plot.

        """
        if self._attractor is None:
            raise ValueError("An attractor must be set up first, please use "
                             "set_attractor(Major_Body).")
        else:
            self._redraw_attractor(trajectory.norm().min() * 0.15)  # Arbitrary threshold

        self._plot_trajectory(trajectory, str(label), color, False)

    def _prepare_plot(self, **layout_kwargs):
        # If there are no orbits, draw only the attractor
        if not self._data:
            self._redraw_attractor()

        if layout_kwargs:
            self._layout.update(layout_kwargs)

    def show(self, **layout_kwargs):
        self._prepare_plot(**layout_kwargs)

        iplot(self.figure)

    def savefig(self, filename, **layout_kwargs):
        self._prepare_plot(**layout_kwargs)

        basename, ext = os.path.splitext(filename)
        export(
            self.figure,
            image=ext[1:], image_filename=basename,
            show_link=False,  # Boilerplate
        )


def plot_solar_system(outer=True, epoch=None):
    """
    Plots the whole solar system in one single call.

    Parameters
    ------------

    outer : To print only the inner solar system , pass outer=False else it will print the whole solar system
    epoch: Epoch Value of the plot. By default the value is None.
    """
    orbits = [Earth, Mars, Mercury, Venus, Jupiter, Saturn, Uranus, Neptune]
    op = OrbitPlotter()
    for orbit in orbits:
        if outer:
            orb = Orbit.from_body_ephem(orbit, epoch)
            op.plot(orb, label=str(orbit))
        else:
            if orbit == Jupiter:
                break
            orb = Orbit.from_body_ephem(orbit, epoch)
            op.plot(orb, label=str(orbit))
    return op
