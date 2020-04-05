import numpy as np
from astropy import units as u
from astropy.coordinates import CartesianRepresentation
from matplotlib import patches as mpl_patches, pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, to_rgba

from poliastro.plotting.util import generate_label

from ._base import BaseOrbitPlotter, Mixin2D, Trajectory


def _segments_from_arrays(x, y):
    # Copied pasted from
    # https://matplotlib.org/3.1.0/gallery/lines_bars_and_markers/multicolored_line.html
    # because this API is impossible to understand
    points = np.column_stack([x.to(u.km).value, y.to(u.km).value])[:, None, :]
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments


class StaticOrbitPlotter(BaseOrbitPlotter, Mixin2D):
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
        super().__init__(num_points)

        self._ax = ax
        if not self._ax:
            if dark:
                with plt.style.context("dark_background"):
                    _, self._ax = plt.subplots(figsize=(6, 6))
            else:
                _, self._ax = plt.subplots(figsize=(6, 6))

        self._frame = None

    def _redraw(self):
        for artist in self._ax.lines + self._ax.collections:
            artist.remove()

        for positions, state, label, colors in self._trajectories:
            self._plot(positions, state, label, colors)

        self._ax.relim()
        self._ax.autoscale()

    def _clear_attractor(self):
        for attractor in self._ax.findobj(match=mpl_patches.Circle):
            attractor.remove()

    def _get_colors(self, color, trail):
        if color is None:
            # HACK: https://stackoverflow.com/a/13831816/554319
            color = next(self._ax._get_lines.prop_cycler)["color"]

        if trail:
            colors = [color, to_rgba(color, 0)]
        else:
            colors = [color]

        return colors

    def _draw_point(self, radius, color, name, center=None):
        x_center, y_center = self._project(
            center[None]
        )  # Indexing trick to add one extra dimension

        (l,) = self._ax.plot(
            x_center.to(u.km).value, y_center.to(u.km).value, "o", mew=0, color=color
        )

        return l

    def _draw_sphere(self, radius, color, name, center=[0, 0, 0] * u.km):
        x_center, y_center = self._project(
            center[None]
        )  # Indexing trick to add one extra dimension

        self._ax.add_patch(
            mpl_patches.Circle(
                (x_center.to(u.km).value, y_center.to(u.km).value),
                radius.to(u.km).value,
                lw=0,
                color=color,
            )
        )

    def _plot_trajectory(self, positions, label, colors, dashed):
        if dashed:
            linestyle = "dashed"
        else:
            linestyle = "solid"

        rr = positions.represent_as(CartesianRepresentation).xyz.transpose()
        x, y = self._project(rr)

        if len(colors) > 1:
            segments = _segments_from_arrays(x, y)
            cmap = LinearSegmentedColormap.from_list(
                f"{colors[0]}_to_alpha", colors  # Useless name
            )
            lc = LineCollection(segments, linestyles=linestyle, cmap=cmap)
            lc.set_array(np.linspace(1, 0, len(x)))

            self._ax.add_collection(lc)
            lines = [lc]

        else:
            lines = self._ax.plot(
                x.to(u.km).value, y.to(u.km).value, linestyle=linestyle, color=colors[0]
            )

        self._ax.set_xlabel("$x$ (km)")
        self._ax.set_ylabel("$y$ (km)")
        self._ax.set_aspect(1)

        return lines

    def _plot_r(self, state, label, colors):
        # TODO: Compute radius?
        return self._draw_point(None, colors[0], label, center=state)

    def _plot(self, positions, state, label, colors):
        trace_trajectory, trace_r = super()._plot(positions, state, label, colors)

        if label:
            if not self._ax.get_legend():
                size = self._ax.figure.get_size_inches() + [8, 0]
                self._ax.figure.set_size_inches(size)

            lines = trace_trajectory[:]
            if trace_r is not None:
                lines.append(trace_r)

            # This will apply the label to either the point or the osculating
            # orbit depending on the last plotted line
            # NOTE: What about generating both labels,
            # indicating that one is the osculating orbit?
            lines[-1].set_label(label)
            self._ax.legend(
                loc="upper left", bbox_to_anchor=(1.05, 1.015), title="Names and epochs"
            )

        return trace_trajectory, trace_r

    def plot_trajectory(self, positions, *, label=None, color=None, trail=False):
        """Plots a precomputed trajectory.

        An attractor must be set first.

        Parameters
        ----------
        positions : ~astropy.coordinates.CartesianRepresentation
            Trajectory to plot.
        label : string, optional
            Label of the trajectory.
        color : string, optional
            Color of the trajectory.
        trail : bool, optional
            Fade the orbit trail, default to False.

        """
        if self._attractor is None:
            raise ValueError(
                "An attractor and a frame must be set up first, please use "
                "set_attractor(Major_Body) or plot(orbit)"
            )
        if self._frame is None:
            raise ValueError(
                "A frame must be set up first, please use "
                "set_orbit_frame(orbit) or plot(orbit)"
            )

        self._redraw_attractor()

        colors = self._get_colors(color, trail)

        # NOTE: We do not call self._plot here to control the labels
        lines = self._plot_trajectory(positions, label, colors, True)

        if label:
            lines[0].set_label(label)
            self._ax.legend(
                loc="upper left", bbox_to_anchor=(1.05, 1.015), title="Names and epochs"
            )

        self._trajectories.append(Trajectory(positions, None, label, colors))

        return lines

    def plot(self, orbit, *, label=None, color=None, trail=False):
        """Plots state and osculating orbit in their plane.

        Parameters
        ----------
        orbit : ~poliastro.twobody.orbit.Orbit
            Orbit to plot.
        label : string, optional
            Label of the orbit.
        color : string, optional
            Color of the line and the position.
        trail : bool, optional
            Fade the orbit trail, default to False.

        """
        if not self._frame:
            self._set_frame(*orbit.pqw())

        colors = self._get_colors(color, trail)

        self.set_attractor(orbit.attractor)
        self._set_plane(orbit.plane, fail_if_set=False)

        if label:
            label = generate_label(orbit.epoch, label)

        positions = orbit.change_plane(self._plane).sample(self._num_points)

        trace_trajectory, trace_r = self._plot(positions, orbit.r, label, colors)

        self._trajectories.append(Trajectory(positions, orbit.r, label, colors))

        lines = trace_trajectory[:]
        if trace_r is not None:
            lines.append(trace_r)

        return lines
