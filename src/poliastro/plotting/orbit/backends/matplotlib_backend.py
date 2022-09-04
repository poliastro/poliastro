"""A module implementing orbit plotter backends based on Matplotlib."""

from astropy import units as u
from matplotlib import patches as mpl_patches, pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, to_rgba

from poliastro.plotting.orbit.backends.base_backend import _OrbitPlotterBackend

def _segments_from_arrays(x, y):
    # Copied pasted from
    # https://matplotlib.org/3.1.0/gallery/lines_bars_and_markers/multicolored_line.html
    # because this API is impossible to understand
    points = np.column_stack([x.to_value(u.km), y.to_value(u.km)])[:, None, :]
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments


class OrbitPlotterBackendMatplotlib2D(_OrbitPlotterBackend):
    """An orbit plotter backend class based on Matplotlib."""

    def __init__(self, ax=None):
        """Initializes a backend instance.

        Parameters
        ----------
        scene : object
            An instance representing the canvas or scene.

        """
        if not ax:
            _, ax = plt.subplots(figsize=(6, 6))
        super().__init__(ax, name="OrbitPlotterBackendMatplotlib2D")

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

        """
        if color is None:
            # HACK: https://stackoverflow.com/a/13831816/554319
            color = next(self.ax._get_lines.prop_cycler)["color"]

        colors = [color, to_rgba(color, 0)] if trail else [color]
        return colors

    def draw_marker(self, position, marker_symbol, *, color=None, size=None):
        """Draws a marker into the scene.

        Parameters
        ----------
        position : list[float, float]
            A list containing the x and y coordinates of the point.
        marker_symbol : str
            The marker symbol to be used when drawing the point.
        color : str, optional
            A string representing the hexadecimal color for the point.
        size : float, optional
            The size of the marker for the dot.

        Returns
        -------
        ~matplotlib.lines.Line2D
            An object representing the trace of the marker in the scene.

        """
        return self.ax.plot(
            position[0].to_value(u.km),
            position[1].to_value(u.km),
            marker=marker_symbol,
            markersize=size,
            color=color,
        )


    def draw_position(self, position, *, color=None, size=None):
        """Draws the position of a body in the scene.

        Parameters
        ----------
        position : list[float, float, float]
            A list containing the x, y and z coordinates of the point.
        color : str, optional
            A string representing the hexadecimal color for the marker.
        size : float, optional
            The size of the marker.

        Returns
        -------
        ~matplotlib.lines.Line2D
            An object representing the trace of the coordinates in the scene.

        """
        return self.draw_marker(position, "o", color=color[0], size=size)


    def draw_impulse(self, position, *, color=None, size=None):
        """Draws an impulse into the scene.

        Parameters
        ----------
        position : list[float, float]
            A list containing the x and y coordinates of the impulse location.
        color : str, optional
            A string representing the hexadecimal color for the impulse marker.
        size : float, optional
            The size of the marker for the impulse.

        Returns
        -------
        ~matplotlib.lines.Line2D
            An object representing the trace of the impulse in the scene.

        """
        return self.draw_marker(position, marker_symbol="x", color=color, size=size)

    def draw_sphere(self, position, *, color=None, size=None):
        """Draws an sphere into the scene.

        Parameters
        ----------
        position : list[float, float]
            A list containing the x and y coordinates of the sphere location.
        color : str, optional
            A string representing the hexadecimal color for the sphere.
        size : float, optional
            The radius of the sphere.

        Returns
        -------
        ~matplotlib.patches.Patch
            An object representing the trace of the sphere in the scene.

        """
        return self.ax.add_patch(
            mpl_patches.Circle(
                (position[0].to_value(u.km), position[1].to_value(u.km)),
                size.to_value(u.km),
                lw=0,
                color=color,
            )
        )

    def undraw_attractor(self):
        """Removes the attractor from the scene."""
        for attractor in self.ax.findobj(match=mpl_patches.Circle):
            attractor.remove()

    def draw_coordinates(self, coordinates, *, color=None, label=None, size=None):
        """Draws desired coordinates into the scene.

        Parameters
        ----------
        position : list[list[float, float], ...]
            A set of lists containing the x and y coordinates of the sphere location.
        color : str, optional
            A string representing the hexadecimal color for the coordinates.
        label : str, optional
            Label of the orbit, default to the name of the body.
        size : float, optional
            The size of the marker for drawing the coordinates.

        Returns
        -------
        trace_coordinates : ~matplotlib.lines.Line2D
            An object representing the trace of the coordinates in the scene.

        """
        # Generate the colors if coordinates trail is required
        if len(color) > 1:
            segments = _segments_from_arrays(x, y)
            cmap = LinearSegmentedColormap.from_list(
                f"{colors[0]}_to_alpha", color  # Useless name
            )
            lc = LineCollection(segments, linestyles=linestyle, cmap=cmap)
            lc.set_array(np.linspace(1, 0, len(x)))

            self.ax.add_collection(lc)
            lines = [lc]

        else:
            # Plot the coordinates in the scene
            trace_coordinates = self.ax.plot(
                coordinates[0].to_value(u.km),
                coordinates[1].to_value(u.km),
                color=color[0],
                label=label,
            )

        # Update the axes labels and impose a 1:1 aspect ratio between them
        self.ax.set_xlabel("$x$ (km)")
        self.ax.set_ylabel("$y$ (km)")
        self.ax.set_aspect(1)

        return trace_coordinates

    def show(self):
        """Displays the scene."""
        plt.show()
