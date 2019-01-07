from typing import List

import matplotlib as mpl
import numpy as np
from astropy import units as u
from astropy.coordinates import CartesianRepresentation
from matplotlib import pyplot as plt

from poliastro.plotting.util import BODY_COLORS, generate_label
from poliastro.twobody.propagation import mean_motion
from poliastro.util import norm

from ._base import Trajectory


class StaticOrbitPlotter:
    """StaticOrbitPlotter class.

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
                with plt.style.context("dark_background"):
                    _, self.ax = plt.subplots(figsize=(6, 6))
            else:
                _, self.ax = plt.subplots(figsize=(6, 6))
        self.num_points = num_points
        self._frame = None
        self._attractor = None
        self._attractor_radius = np.inf * u.km
        self._trajectories = []  # type: List[Trajectory]

    @property
    def trajectories(self):
        return self._trajectories

    def set_frame(self, p_vec, q_vec, w_vec):
        """Sets perifocal frame.

        Raises
        ------
        ValueError
            If the vectors are not a set of mutually orthogonal unit vectors.
        """
        if not np.allclose([norm(v) for v in (p_vec, q_vec, w_vec)], 1):
            raise ValueError("Vectors must be unit.")
        elif not np.allclose([p_vec.dot(q_vec), q_vec.dot(w_vec), w_vec.dot(p_vec)], 0):
            raise ValueError("Vectors must be mutually orthogonal.")
        else:
            self._frame = p_vec, q_vec, w_vec

        if self._trajectories:
            self._redraw()

    def _redraw(self):
        for artist in self.ax.lines + self.ax.collections:
            artist.remove()

        for trajectory, state, label, color in self._trajectories:
            self._plot(trajectory, state, label, color)

        self.ax.relim()
        self.ax.autoscale()

    def _plot_trajectory(self, trajectory, color=None):
        rr = trajectory.represent_as(CartesianRepresentation).xyz.transpose()
        x, y = self._project(rr)
        lines = self.ax.plot(x.to(u.km).value, y.to(u.km).value, "--", color=color)

        return lines

    def plot_trajectory(self, trajectory, *, label=None, color=None):
        """Plots a precomputed trajectory.

        Parameters
        ----------
        trajectory : ~astropy.coordinates.BaseRepresentation, ~astropy.coordinates.BaseCoordinateFrame
            Trajectory to plot.
        label : str, optional
            Label.
        color : str, optional
            Color string.

        """
        lines = self._plot_trajectory(trajectory, color)

        if label:
            lines[0].set_label(label)
            self.ax.legend(
                loc="upper left", bbox_to_anchor=(1.05, 1.015), title="Names and epochs"
            )

        self._trajectories.append(
            Trajectory(trajectory, None, label, lines[0].get_color())
        )

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
            raise NotImplementedError(
                "Attractor has already been set to {}.".format(self._attractor.name)
            )

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

        self.ax.add_patch(
            mpl.patches.Circle((0, 0), self._attractor_radius.value, lw=0, color=color)
        )

    def _plot(self, trajectory, state=None, label=None, color=None):
        lines = self._plot_trajectory(trajectory, color)

        if state is not None:
            x0, y0 = self._project(state[None])

            # Plot current position
            l, = self.ax.plot(
                x0.to(u.km).value,
                y0.to(u.km).value,
                "o",
                mew=0,
                color=lines[0].get_color(),
            )
            lines.append(l)

        if label:
            if not self.ax.get_legend():
                size = self.ax.figure.get_size_inches() + [8, 0]
                self.ax.figure.set_size_inches(size)

            # This will apply the label to either the point or the osculating
            # orbit depending on the last plotted line
            lines[-1].set_label(label)
            self.ax.legend(
                loc="upper left", bbox_to_anchor=(1.05, 1.015), title="Names and epochs"
            )

        self.ax.set_xlabel("$x$ (km)")
        self.ax.set_ylabel("$y$ (km)")
        self.ax.set_aspect(1)

        return lines

    def plot(self, orbit, label=None, color=None, method=mean_motion):
        """Plots state and osculating orbit in their plane.
        """
        if not self._frame:
            self.set_frame(*orbit.pqw())

        self.set_attractor(orbit.attractor)
        self._redraw_attractor(orbit.r_p * 0.15)  # Arbitrary threshold
        positions = orbit.sample(self.num_points, method)

        if label:
            label = generate_label(orbit, label)

        lines = self._plot(positions, orbit.r, label, color)

        self._trajectories.append(
            Trajectory(positions, orbit.r, label, lines[0].get_color())
        )

        return lines
