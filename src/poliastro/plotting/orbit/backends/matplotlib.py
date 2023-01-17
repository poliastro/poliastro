"""A module implementing orbit plotter backends based on Matplotlib."""

import numpy as np
from matplotlib import patches as mpl_patches, pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, to_rgba

from poliastro.plotting.orbit.backends._base import OrbitPlotterBackend


def _segments_from_arrays(x, y):
    # Copied pasted from
    # https://matplotlib.org/3.1.0/gallery/lines_bars_and_markers/multicolored_line.html
    # because this API is impossible to understand
    points = np.column_stack([x, y])[:, None, :]
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments


class Matplotlib2D(OrbitPlotterBackend):
    """An orbit plotter backend class based on Matplotlib."""

    def __init__(self, ax=None, use_dark_theme=False):
        """Initializes a backend instance.

        Parameters
        ----------
        ax: ~matplotlib.axes.Axes
            An :ref:`~matplotlib.axes.Axes` instance representing the axes of the figure.
        use_dark_theme : bool, optional
            If ``True``, uses dark theme. If ``False``, uses light theme.
            Default to ``False``.

        """
        if ax is None:
            if use_dark_theme is True:
                with plt.style.context("dark_background"):
                    _, ax = plt.subplots(figsize=(6, 6))
            else:
                _, ax = plt.subplots(figsize=(6, 6))
        super().__init__(ax, self.__class__.__name__)

    @property
    def ax(self):
        """The matplotlib axes were the scene is rendered.

        Returns
        -------
        ax: ~matplotlib.axes.Axes
            The matplotlib Axes representing the scene.

        """
        return self.scene

    def _get_colors(self, color, trail):
        """Return the required list of colors if orbit trail is desired.

        Parameters
        ----------
        color : str
            A string representing the hexadecimal color for the point.
        trail : bool
            ``True`` if orbit trail is desired, ``False`` if not desired.

        Returns
        -------
        list[str]
            A list of strings representing hexadecimal colors.

        """
        if color is None:
            # HACK: https://stackoverflow.com/a/13831816/554319
            color = next(self.ax._get_lines.prop_cycler)["color"]

        colors = [color, to_rgba(color, 0)] if trail else [color]
        return colors

    def draw_marker(self, position, *, color, label, marker_symbol, size):
        """Draw a marker into the scene.

        Parameters
        ----------
        position : list[float, float]
            A list containing the x and y coordinates of the point.
        color : str
            A string representing the hexadecimal color for the point.
        label : str
            The name to be used in the legend for the marker.
        marker_symbol : str
            The marker symbol to be used when drawing the point.
        size : float
            Desired size for the marker.

        Returns
        -------
        ~matplotlib.lines.Line2D
            An object representing the trace of the marker in the scene.

        """
        (marker_trace,) = self.ax.plot(
            position[0],
            position[1],
            color=color,
            marker=marker_symbol,
            markersize=size,
            markeredgecolor=color,
            label=label,
            linestyle="None",
        )
        return marker_trace

    def draw_position(self, position, *, color, label, size):
        """Draw the position of a body in the scene.

        Parameters
        ----------
        position : list[float, float, float]
            A list containing the x, y and z coordinates of the point.
        color : str
            A string representing the hexadecimal color for the marker.
        label : str
            The name to be used to identify the position in the legend of the
            figure.
        size : float
            The size of the marker.

        Returns
        -------
        ~matplotlib.lines.Line2D
            An object representing the trace of the coordinates in the scene.

        """
        return self.draw_marker(
            position, color=color, marker_symbol="o", label=label, size=None
        )

    def draw_impulse(self, position, *, color, label, size):
        """Draw an impulse into the scene.

        Parameters
        ----------
        position : list[float, float]
            A list containing the x and y coordinates of the impulse location.
        color : str
            A string representing the hexadecimal color for the impulse marker.
        label : str
            The name to be used to identify the position in the legend of the
            figure.
        size : float
            The size of the marker for the impulse.

        Returns
        -------
        ~matplotlib.lines.Line2D
            An object representing the trace of the impulse in the scene.

        """
        return self.draw_marker(
            position, color=color, label=label, marker_symbol="x", size=size
        )

    def draw_sphere(self, position, *, color, label, radius):
        """Draw an sphere into the scene.

        Parameters
        ----------
        position : list[float, float]
            A list containing the x and y coordinates of the sphere location.
        color : str
            A string representing the hexadecimal color for the sphere.
        label : str
            The name shown in the legend of the figure to identify the sphere.
        radius : float
            The radius of the sphere.

        Returns
        -------
        ~matplotlib.patches.Patch
            An object representing the trace of the sphere in the scene.

        """
        return self.ax.add_patch(
            mpl_patches.Circle(
                (
                    position[0],
                    position[1],
                ),
                radius,
                color=color,
                linewidth=0,
                label=label,
            )
        )

    def undraw_attractor(self):
        """Remove the attractor from the scene."""
        for attractor in self.ax.findobj(match=mpl_patches.Circle):
            attractor.remove()

    def draw_axes_labels_with_length_scale_units(self, length_scale_units):
        """Draw the desired label into the specified axis.

        Parameters
        ----------
        lenght_scale_units : ~astropy.units.Unit
            Desired units of lenght used for representing distances.

        """
        self.ax.set_xlabel(f"$x$ ({length_scale_units.name})")
        self.ax.set_ylabel(f"$y$ ({length_scale_units.name})")

    def draw_coordinates(self, coordinates, *, colors, dashed, label):
        """Draw desired coordinates into the scene.

        Parameters
        ----------
        coordinates : list[list[float, float, float]]
            A set of lists containing the x, y and z coordinates.
        colors : list[str]
            A list of string representing the hexadecimal color for the coordinates.
        dashed : bool
            Whether to use a dashed or solid line style for the coordiantes.
        label : str
            The name to be used to identify the coordinates in the legend of the
            figure.

        Returns
        -------
        trace_coordinates : ~matplotlib.lines.Line2D
            An object representing the trace of the coordinates in the scene.

        """
        # Select the desired linestyle for the line representing the coordinates
        linestyle = "dashed" if dashed else "solid"

        # Unpack the x and y coordinates
        x, y, _ = coordinates

        # Generate the colors if coordinates trail is required
        if len(colors) > 1:
            segments = _segments_from_arrays(x, y)
            cmap = LinearSegmentedColormap.from_list(
                f"{colors[0]}_to_alpha", colors  # Useless name
            )
            lc = LineCollection(
                segments, linestyles=linestyle, cmap=cmap, label=label
            )
            lc.set_array(np.linspace(1, 0, len(x)))

            self.ax.add_collection(lc)
            lines_coordinates = [lc]

        else:
            # Plot the coordinates in the scene
            (lines_coordinates,) = self.ax.plot(
                x,
                y,
                color=colors[0],
                label=label,
                linestyle=linestyle,
            )

        return lines_coordinates

    def generate_labels(self, label, has_coordinates, has_position):
        """Generate the labels for coordinates and position.

        Parameters
        ----------
        label : str
            A string representing the label.
        has_coordinates : boolean
            Whether the object has coordinates to plot or not.
        has_position : boolean
            Whether the object has a position to plot or not.

        Returns
        -------
        tuple
            A tuple containing the coordinates and position labels.

        """
        return (None, label) if has_position else (label, None)

    def update_legend(self):
        """Update the legend of the scene."""
        # Enable the legend (if required)
        if not self.ax.get_legend():
            size = self.ax.figure.get_size_inches() + [8, 0]
            self.ax.figure.set_size_inches(size)

        self.ax.legend(
            loc="upper left",
            bbox_to_anchor=(1.05, 1.015),
            title="Names and epochs",
            numpoints=1,
        )

    def resize_limits(self):
        """Resize the limits of the scene."""
        self.ax.relim()
        self.ax.autoscale()
        self.ax.set_aspect(1)

    def show(self):
        """Display the scene."""
        plt.show()
