import numpy as np
from astropy import units as u
from matplotlib import patches as mpl_patches, pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, to_rgba

from ._base import BaseOrbitPlotter, Mixin2D


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

    def __init__(self, ax=None, num_points=150, dark=False, *, plane=None):
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
        super().__init__(num_points=num_points, plane=plane)

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

        super()._redraw()

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

    def _plot_coordinates(self, coordinates, label, colors, dashed):
        if self._frame is None:
            raise ValueError(
                "A frame must be set up first, please use "
                "set_orbit_frame(orbit) or plot(orbit)"
            )

        if dashed:
            linestyle = "dashed"
        else:
            linestyle = "solid"

        rr = coordinates.xyz.transpose()
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

    def _plot_position(self, position, label, colors):
        # TODO: Compute radius?
        return self._draw_point(None, colors[0], label, center=position)

    def _set_legend(self, label, line_coordinates, line_position=None):
        if not self._ax.get_legend():
            size = self._ax.figure.get_size_inches() + [8, 0]
            self._ax.figure.set_size_inches(size)

        # This will apply the label to either the point or the osculating
        # orbit depending on the last plotted line
        if line_position is not None:
            line_position.set_label(label)
        else:
            line_coordinates[0].set_label(label)

        self._ax.legend(
            loc="upper left",
            bbox_to_anchor=(1.05, 1.015),
            title="Names and epochs",
            numpoints=1,
        )

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

        lines = self._plot_trajectory(
            coordinates, label=label, color=color, trail=trail
        )

        if label:
            self._set_legend(label, *lines)

        return lines

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
        if not self._frame:
            self.set_orbit_frame(orbit)

        lines = self._plot(orbit, label=label, color=color, trail=trail)
        lines = lines[0] + [lines[1]]

        # Set legend using label from last added trajectory
        self._set_legend(self._trajectories[-1].label, *lines)

        return lines

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
        if self._frame is None:
            self.set_body_frame(body, epoch)

        lines = self._plot_body_orbit(
            body, epoch, label=label, color=color, trail=trail
        )

        # Set legend using label from last added trajectory
        self._set_legend(self._trajectories[-1].label, *lines)

        return lines

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

        lines = self._plot_ephem(ephem, epoch, label=label, color=color, trail=trail)

        self._set_legend(self._trajectories[-1].label, *lines)

        return lines
