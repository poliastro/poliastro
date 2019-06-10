from typing import List

import matplotlib as mpl
import numpy as np
from astropy import units as u
from astropy.coordinates import CartesianRepresentation
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, to_rgba

from poliastro.plotting.util import BODY_COLORS, generate_label

from ._base import Trajectory


def _segments_from_arrays(x, y):
    # Copied pasted from
    # https://matplotlib.org/3.1.0/gallery/lines_bars_and_markers/multicolored_line.html
    # because this API is impossible to understand
    points = np.column_stack([x.to(u.km).value, y.to(u.km).value])[:, None, :]
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments


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

    def set_frame(self, frame):
        """Sets perifocal frame.

        """
        self._frame = frame

        if self._trajectories:
            self._redraw()

    def _get_colors(self, color, trail):
        if trail and color is None:
            # HACK: https://stackoverflow.com/a/13831816/554319
            color = next(self.ax._get_lines.prop_cycler)["color"]

        if trail:
            colors = [color, to_rgba(color, 0)]
        else:
            colors = [color]

        return colors

    def _redraw(self):
        for artist in self.ax.lines + self.ax.collections:
            artist.remove()

        for trajectory, state, label, colors in self._trajectories:
            self._plot(trajectory, state, label, colors)

        self.ax.relim()
        self.ax.autoscale()

    def _plot_trajectory(self, trajectory, colors=None, linestyle="dashed"):
        x, y = self._project(trajectory)

        if len(colors) > 1:
            segments = _segments_from_arrays(x, y)
            cmap = LinearSegmentedColormap.from_list(
                f"{colors[0]}_to_alpha", colors  # Useless name
            )
            lc = LineCollection(segments, linestyles=linestyle, cmap=cmap)
            lc.set_array(np.linspace(1, 0, len(x)))

            self.ax.add_collection(lc)
            lines = [lc]

        else:
            lines = self.ax.plot(
                x.to(u.km).value, y.to(u.km).value, linestyle=linestyle, color=colors[0]
            )
            colors = [lines[0].get_color()]

        return lines, colors

    def plot_trajectory(self, trajectory, *, label=None, color=None, trail=False):
        """Plots a precomputed trajectory.

        Parameters
        ----------
        trajectory : ~astropy.coordinates.BaseRepresentation, ~astropy.coordinates.BaseCoordinateFrame
            Trajectory to plot.
        label : str, optional
            Label.
        color : str, optional
            Color string.
        trail: bool, optional
            Plots the Orbit's trail

        """
        if self._attractor is None or self._frame is None:
            raise ValueError(
                "An attractor and a frame must be set up first, please use "
                "set_attractor(Major_Body) and set_frame(orbit.get_perifocal_frame()) "
                "or plot(orbit)."
            )

        # Strip velocities from trajectory if present
        # HACK: Is there a better way?
        trajectory = trajectory.copy()
        trajectory.data.differentials.pop("s", None)

        self._redraw_attractor(
            trajectory.represent_as(CartesianRepresentation).norm().min() * 0.15
        )  # Arbitrary threshold

        colors = self._get_colors(color, trail)
        lines, colors = self._plot_trajectory(trajectory, colors)

        if label:
            lines[0].set_label(label)
            self.ax.legend(
                loc="upper left", bbox_to_anchor=(1.05, 1.015), title="Names and epochs"
            )

        self._trajectories.append(Trajectory(trajectory, None, label, colors))

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
                f"Attractor has already been set to {self._attractor.name}."
            )

    def _project(self, rr):
        rr_proj = rr.transform_to(self._frame).represent_as(CartesianRepresentation)
        return rr_proj.x, rr_proj.y

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

    def _plot(self, trajectory, state=None, label=None, colors=None):
        lines, colors = self._plot_trajectory(trajectory, colors)

        if state is not None:
            x0, y0 = self._project(state[None])

            # Plot current position
            (l,) = self.ax.plot(
                x0.to(u.km).value, y0.to(u.km).value, "o", mew=0, color=colors[0]
            )
            lines.append(l)

        if label:
            if not self.ax.get_legend():
                size = self.ax.figure.get_size_inches() + [8, 0]
                self.ax.figure.set_size_inches(size)

            # This will apply the label to either the point or the osculating
            # orbit depending on the last plotted line
            # NOTE: What about generating both labels,
            # indicating that one is the osculating orbit?
            lines[-1].set_label(label)
            self.ax.legend(
                loc="upper left", bbox_to_anchor=(1.05, 1.015), title="Names and epochs"
            )

        self.ax.set_xlabel("$x$ (km)")
        self.ax.set_ylabel("$y$ (km)")
        self.ax.set_aspect(1)

        return lines, colors

    def plot(self, orbit, label=None, color=None, trail=False):
        """Plots state and osculating orbit in their plane.

        """
        if not self._frame:
            self.set_frame(orbit.get_perifocal_frame())

        self.set_attractor(orbit.attractor)
        self._redraw_attractor(orbit.r_p * 0.15)  # Arbitrary threshold

        positions = orbit.sample(self.num_points)

        if label:
            label = generate_label(orbit, label)

        # Use of a protected method instead of frame.realize_frame
        # because the latter does not let the user choose the representation type
        # in one line despite its parameter names, see
        # https://github.com/astropy/astropy/issues/7784
        r_framed = orbit.frame._replicate(orbit.r, representation_type="cartesian")

        colors = self._get_colors(color, trail)
        lines, colors = self._plot(positions, r_framed, label, colors)

        self._trajectories.append(Trajectory(positions, r_framed, label, colors))
        return lines
