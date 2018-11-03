""" Plotting utilities.

"""
import numpy as np

from typing import List, Tuple  # flake8: noqa

import matplotlib as mpl
import matplotlib.pyplot as plt

from astropy import units as u
from astropy.coordinates import CartesianRepresentation

from poliastro.util import norm
from poliastro.twobody.propagation import mean_motion


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


class OrbitPlotter(object):
    """OrbitPlotter class.

    This class holds the perifocal plane of the first
    :py:class:`~poliastro.twobody.orbit.Orbit` plotted in it using
    :py:meth:`plot`, so all following
    plots will be projected on that plane. Alternatively, you can call
    :py:meth:`set_frame` to set the frame before plotting.

    """

    def __init__(self, ax=None, num_points=150, dark=False):
        """Constructor.

        Parameters
        ----------
        ax : ~matplotlib.axes.Axes
            Axes in which to plot. If not given, new ones will be created.
        num_points : int, optional
            Number of points to use in plots, default to 150.
        dark : bool, optional
            If set as True, plots the orbit in Dark mode.
        """
        self.ax = ax
        if not self.ax:
            if dark:
                with plt.style.context('dark_background'):
                    _, self.ax = plt.subplots(figsize=(6, 6))
            else:
                _, self.ax = plt.subplots(figsize=(6, 6))
        self.num_points = num_points
        self._frame = None
        self._attractor = None
        self._attractor_radius = np.inf * u.km
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
        self._attractor = None
        for orbit, label in self._orbits:
            self.plot(orbit, label)
        self.ax.relim()
        self.ax.autoscale()

    def plot_trajectory(self, trajectory, *, label=None, color=None):
        """Plots a precomputed trajectory.

        Parameters
        ----------
        trajectory : ~astropy.coordinates.BaseRepresentation, ~astropy.coordinates.BaseCoordinateFrame
            Trajectory to plot.

        """
        lines = []
        rr = trajectory.represent_as(CartesianRepresentation).xyz.transpose()
        x, y = self._project(rr)
        a, = self.ax.plot(x.to(u.km).value, y.to(u.km).value, '--', color=color, label=label)
        lines.append(a)
        if label:
            a.set_label(label)
            self.ax.legend(bbox_to_anchor=(1.05, 1), title="Names and epochs")

        return lines

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

    def _project(self, rr):
        rr_proj = rr - rr.dot(self._frame[2])[:, None] * self._frame[2]
        x = rr_proj.dot(self._frame[0])
        y = rr_proj.dot(self._frame[1])
        return x, y

    def _redraw_attractor(self, min_radius=0 * u.km):
        radius = max(self._attractor.R.to(u.km), min_radius.to(u.km))
        color = BODY_COLORS.get(self._attractor.name, "#999999")

        for attractor in self.ax.findobj(match=mpl.patches.Circle):
            attractor.remove()

        if radius < self._attractor_radius:
            self._attractor_radius = radius

        self.ax.add_patch(mpl.patches.Circle((0, 0), self._attractor_radius.value, lw=0, color=color))

    def plot(self, orbit, label=None, color=None, method=mean_motion):
        """Plots state and osculating orbit in their plane.
        """
        if not self._frame:
            self.set_frame(*orbit.pqw())

        if (orbit, label) not in self._orbits:
            self._orbits.append((orbit, label))

        self.set_attractor(orbit.attractor)
        self._redraw_attractor(orbit.r_p * 0.15)  # Arbitrary Threshhold
        positions = orbit.sample(self.num_points, method)

        x0, y0 = self._project(orbit.r[None])
        # Plot current position
        l, = self.ax.plot(x0.to(u.km).value, y0.to(u.km).value,
                          'o', mew=0, color=color)

        lines = self.plot_trajectory(trajectory=positions, color=l.get_color())
        lines.append(l)

        if label:
            # This will apply the label to either the point or the osculating
            # orbit depending on the last plotted line, as they share variable
            if not self.ax.get_legend():
                size = self.ax.figure.get_size_inches() + [8, 0]
                self.ax.figure.set_size_inches(size)
            label = _generate_label(orbit, label)
            l.set_label(label)
            self.ax.legend(bbox_to_anchor=(1.05, 1), title="Names and epochs")

        self.ax.set_xlabel("$x$ (km)")
        self.ax.set_ylabel("$y$ (km)")
        self.ax.set_aspect(1)

        return lines


def plot(state, label=None, color=None):
    """Plots an :py:class:`~poliastro.twobody.orbit.Orbit` in 2D.

    For more advanced tuning, use the :py:class:`OrbitPlotter` class.

    """
    op = OrbitPlotter()
    return op.plot(state, label=label, color=color)
