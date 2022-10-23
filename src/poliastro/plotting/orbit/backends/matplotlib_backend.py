"""A module implementing orbit plotter backends based on Matplotlib."""

import numpy as np
from astropy import units as u
from matplotlib import patches as mpl_patches, pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, to_rgba

from poliastro.plotting.orbit.backends.base_backend import _OrbitPlotterBackend


def _segments_from_arrays(x, y):
    # Copied pasted from
    # https://matplotlib.org/3.1.0/gallery/lines_bars_and_markers/multicolored_line.html
    # because this API is impossible to understand
    points = np.column_stack([x, y])[:, None, :]
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments


class OrbitPlotterBackendMatplotlib2D(_OrbitPlotterBackend):
    """An orbit plotter backend class based on Matplotlib."""

    def __init__(self, ax=None, use_dark_theme=None):
        """Initializes a backend instance.

        Parameters
        ----------
        scene : object
            An instance representing the canvas or scene.
        use_dark_theme : bool, optional
            If ``True``, uses dark theme. If ``False``, uses light theme.
            Default to ``False``.

        """
        if not ax:
            if use_dark_theme is True:
                with plt.style.context("dark_background"):
                    _, ax = plt.subplots(figsize=(6, 6))
            else:
                _, ax = plt.subplots(figsize=(6, 6))
        super().__init__(ax, name=self.__class__.__name__)

    @property
    def ax(self):
        """Return the Matplotlib axes were the scene is rendered.

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

    def draw_label(self, label, trace_coordinates, trace_position):
        """Draw the desired label in the figure's legend.

        Parameters
        ----------
        label : str
             A string representing the label name to be drawn in the legend.
        trace_coordinates : object
            An object representing the trace of the coordinates in the scene.
        trace_position : object
            An object representing the trace of the position in the scene.

        """
        if not self.ax.get_legend():
            size = self.ax.figure.get_size_inches() + [8, 0]
            self.ax.figure.set_size_inches(size)

        # This will apply the label to either the point or the osculating
        # orbit depending on the last plotted line
        trace = trace_position or trace_coordinates
        trace[0].set_label(label)

        self.ax.legend(
            loc="upper left",
            bbox_to_anchor=(1.05, 1.015),
            title="Names and epochs",
            numpoints=1,
        )

    def draw_marker(self, position, *, color, marker_symbol, label, size):
        """Draws a marker into the scene.

        Parameters
        ----------
        position : list[float, float]
            A list containing the x and y coordinates of the point.
        color : str, optional
            A string representing the hexadecimal color for the point.
        marker_symbol : str
            The marker symbol to be used when drawing the point.
        label : str
            A string containing tha message for the label.
        size : float, optional
            Desired size for the marker.

        Returns
        -------
        ~matplotlib.lines.Line2D
            An object representing the trace of the marker in the scene.

        """
        return self.ax.plot(
            position[0].to_value(u.km),
            position[1].to_value(u.km),
            color=color,
            marker=marker_symbol,
            markersize=size,
            label=label,
            linestyle="None",
        )

    def draw_position(self, position, *, color, label, size):
        """Draws the position of a body in the scene.

        Parameters
        ----------
        position : list[float, float, float]
            A list containing the x, y and z coordinates of the point.
        color : str, optional
            A string representing the hexadecimal color for the marker.
        size : float, optional
            The size of the marker.
        label : str
            A string containing tha message for the label.

        Returns
        -------
        ~matplotlib.lines.Line2D
            An object representing the trace of the coordinates in the scene.

        """
        return self.draw_marker(
            position, color=color, marker_symbol="o", label=None, size=size
        )

    def draw_impulse(self, position, *, color, label, size):
        """Draws an impulse into the scene.

        Parameters
        ----------
        position : list[float, float]
            A list containing the x and y coordinates of the impulse location.
        color : str, optional
            A string representing the hexadecimal color for the impulse marker.
        label : str
            A string containing tha message for the label.
        size : float, optional
            The size of the marker for the impulse.

        Returns
        -------
        ~matplotlib.lines.Line2D
            An object representing the trace of the impulse in the scene.

        """
        return self.draw_marker(
            position, color=color, label=label, marker_symbol="x", size=size
        )

    def draw_sphere(self, position, *, color, radius):
        """Draws an sphere into the scene.

        Parameters
        ----------
        position : list[float, float]
            A list containing the x and y coordinates of the sphere location.
        color : str, optional
            A string representing the hexadecimal color for the sphere.
        radius : float, optional
            The radius of the sphere.

        Returns
        -------
        ~matplotlib.patches.Patch
            An object representing the trace of the sphere in the scene.

        """
        return self.ax.add_patch(
            mpl_patches.Circle(
                (position[0].to_value(u.km), position[1].to_value(u.km)),
                radius.to_value(u.km),
                color=color,
                linewidth=0,
            )
        )

    def undraw_attractor(self):
        """Removes the attractor from the scene."""
        for attractor in self.ax.findobj(match=mpl_patches.Circle):
            attractor.remove()

    def draw_coordinates(self, coordinates, *, colors, dashed):
        """Draws desired coordinates into the scene.

        Parameters
        ----------
        position : list[list[float, float, float]]
            A set of lists containing the x and y coordinates of the sphere location.
        colors : list[str]
            A list of string representing the hexadecimal color for the coordinates.
        dashed : bool
            Whether to use a dashed or solid line style for the coordiantes.
        label : str
            Label of the orbit, default to the name of the body.

        Returns
        -------
        trace_coordinates : ~matplotlib.lines.Line2D
            An object representing the trace of the coordinates in the scene.

        """
        # Select the desired linestyle for the line representing the coordinates
        linestyle = "dashed" if dashed else "solid"

        # Unpack coordinates
        x, y, z = (coords.to_value(u.km) for coords in coordinates)

        # Generate the colors if coordinates trail is required
        if len(colors) > 1:
            segments = _segments_from_arrays(x, y)
            cmap = LinearSegmentedColormap.from_list(
                f"{colors[0]}_to_alpha", colors  # Useless name
            )
            lc = LineCollection(segments, linestyles=linestyle, cmap=cmap)
            lc.set_array(np.linspace(1, 0, len(x)))

            self.ax.add_collection(lc)
            lines_coordinates = [lc]

        else:
            # Plot the coordinates in the scene
            lines_coordinates = self.ax.plot(
                x,
                y,
                color=colors[0],
                linestyle=linestyle,
            )

        # Update the axes labels and impose a 1:1 aspect ratio between them
        self.ax.set_xlabel("$x$ (km)")
        self.ax.set_ylabel("$y$ (km)")
        self.ax.set_aspect(1)

        return lines_coordinates

    def show(self):
        """Displays the scene."""
        plt.show()
