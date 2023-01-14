"""A module implementing orbit plotter backends based on Plotly."""

from itertools import cycle

import plotly
import plotly.graph_objects as go

from poliastro.plotting.orbit.backends._base import OrbitPlotterBackend
from poliastro.plotting.util import generate_sphere


class BasePlotly(OrbitPlotterBackend):
    """An orbit plotter backend class based on Plotly."""

    def __init__(self, figure, layout):
        """Initializes a backend instance.

        Parameters
        ----------
        figure : ~plotly.graph_objects.Figure
            The plotly ``Figure`` to render the scene.
        layout : ~plotly.graph_objects.Layout
            The plotly ``Layout`` object linked to the figure.

        """
        figure = figure or go.Figure()
        super().__init__(figure, self.__class__.__name__)

        self._layout = layout or go.Layout()
        self.update_layout(self._layout)

        self._color_cycle = cycle(plotly.colors.DEFAULT_PLOTLY_COLORS)

    @property
    def figure(self):
        """Return the Matplotlib axes were the scene is rendered.

        Returns
        -------
        ~plotly.graph_objects.Figure
            The plotly ``Figure`` representing the scene.

        """
        return self.scene

    @property
    def layout(self):
        """Return the layout of the figure.

        Returns
        -------
        ~plotly.graph_objects.Layout
            The plotly ``Layout`` object linked to the figure.

        """
        return self._layout

    def update_layout(self, layout):
        """Update the layout of the figure scene.

        Parameters
        ----------
        layout : ~plotly.graph_objects.Layout
            The new plotly ``Layout`` to be used in the figure.

        """
        self.figure.update_layout(layout)

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
        color = color or next(self._color_cycle)
        return [color]

    def undraw_attractor(self):
        """Removes the attractor from the scene."""
        pass

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
            The name shown in the figure legend to identify the position.

        Returns
        -------
        [~plotly.graph_objects.Surface, ~plotly.graph_objects.Trace]
            An object representing the trace of the coordinates in the scene.

        """
        return self.draw_sphere(
            position, color=color, label=label, radius=size
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
            The name shown in the figure legend to identify the impulse.
        size : float, optional
            The size of the marker for the impulse.

        Returns
        -------
        object
            An object representing the trace of the impulse in the scene.

        """
        return self.draw_marker(
            position, color=color, label=label, marker_symbol="x", size=size
        )

    def update_legend(self):
        """Update the legend of the scene."""
        pass

    def resize_limits(self):
        """Resize the limits of the scene."""
        pass

    def show(self):
        """Displays the scene."""
        self.update_layout(self._layout)
        if not self.figure._in_batch_mode:
            return self.figure.show()

    def generate_labels(self, label, has_coordinates, has_position):
        """Generates the labels for coordinates and position.

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
        return (label, None)


class Plotly2D(BasePlotly):
    """An orbit plotter backend class based on Plotly."""

    def __init__(self, figure, use_dark_theme):
        """Initializes a backend instance.

        Parameters
        ----------
        figure : ~plotly.graph_objects.Figure
            The plotly ``Figure`` to render the scene.
        use_dark_theme : bool, optional
            If ``True``, uses dark theme. If ``False``, uses light theme.
            Default to ``False``.

        """
        # Apply the desired theme
        theme = "plotly_dark" if use_dark_theme is True else "plotly"

        # Declare the layout and attach it to the figure
        layout = go.Layout(
            autosize=True,
            xaxis=dict(constrain="domain"),
            yaxis=dict(scaleanchor="x"),
            template=theme,
        )
        super().__init__(figure, layout)

    def draw_marker(self, position, *, color, label, marker_symbol, size):
        """Draws a marker into the scene.

        Parameters
        ----------
        position : list[float, float]
            A list containing the x and y coordinates of the point.
        color : str, optional
            A string representing the hexadecimal color for the point.
        label : str
            The name shown in the legend of the figure to identify the marker.
        marker_symbol : str
            The marker symbol to be used when drawing the point.
        size : float, optional
            Desired size for the marker.

        Returns
        -------
        object
            An object representing the trace of the marker in the scene.

        """
        marker_style = dict(size=size, color=color, symbol=marker_symbol)
        marker_trace = go.Scatter(
            x=position[0],
            y=position[1],
            marker=marker_style,
            name=label,
            showlegend=False if label is None else True,
        )
        self.figure.add_trace(marker_trace)
        return marker_trace

    def draw_sphere(self, position, *, color, label, radius):
        """Draws an sphere into the scene.

        Parameters
        ----------
        position : list[float, float]
            A list containing the x and y coordinates of the sphere location.
        color : str, optional
            A string representing the hexadecimal color for the sphere.
        label : str
            Unuseful for this routine. See the ``Notes`` section.
        radius : float, optional
            The radius of the sphere.

        Notes
        -----
        Plotting a sphere in a two-dimensional figure in plotly requires a shape
        instead of a trace. Shapes do not accept a label, as the legend does not
        support labels for shapes.

        Returns
        -------
        dict
            A dictionary representing the shape of the sphere.

        """
        shape = dict(
            type="circle",
            xref="x",
            yref="y",
            x0=(position[0] - radius),
            y0=(position[1] - radius),
            x1=(position[0] + radius),
            y1=(position[1] + radius),
            fillcolor=color,
            line=dict(color=color),
            opacity=1,
        )
        self.layout.shapes += (shape,)
        return shape

    def draw_coordinates(self, coordinates, *, colors, dashed, label):
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
            The name shown in the legend for identifying the coordinates.

        Returns
        -------
        trace_coordinates : object
            An object representing the trace of the coordinates in the scene.

        """
        # Select the desired linestyle for the line representing the coordinates
        linestyle = "dash" if dashed else "solid"

        # Unpack coordinates
        x, y, _ = coordinates

        # Plot the coordinates in the scene
        coordinates_trace = go.Scatter(
            x=x,
            y=y,
            line=dict(color=colors[0], width=5, dash=linestyle),
            mode="lines",
            name=label,
            showlegend=False if label is None else True,
        )
        self.figure.add_trace(coordinates_trace)
        return coordinates_trace

    def draw_axes_labels_with_length_scale_units(self, length_scale_units):
        """Draws the desired label into the specified axis.

        Parameters
        ----------
        lenght_scale_units : ~astropy.units.Unit
            Desired units of lenght used for representing distances.

        """
        # HACK: plotly does not show LaTeX symbols and \text. The usage of
        # ASCII labels in plotly figures is used to avoid this issue
        self.figure.update_layout(
            xaxis_title=f"x ({length_scale_units.name})",
            yaxis_title=f"y ({length_scale_units.name})",
        )


class Plotly3D(BasePlotly):
    """An orbit plotter backend class based on Plotly."""

    def __init__(self, figure, use_dark_theme):
        """Initializes a backend instance.

        Parameters
        ----------
        figure : ~plotly.graph_objects.Figure
            The plotly ``Figure`` to render the scene.
        use_dark_theme : bool, optional
            If ``True``, uses dark theme. If ``False``, uses light theme.
            Default to ``False``.

        """
        # Apply the desired theme
        theme = "plotly_dark" if use_dark_theme is True else "plotly"

        # Declare the layout and attach it to the figure
        layout = go.Layout(
            autosize=True,
            scene=dict(
                aspectmode="data",
            ),
            template=theme,
        )
        super().__init__(figure, layout)

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
            The name shown in the legend of the figure to identify the marker.
        size : float, optional
            Desired size for the marker.

        Returns
        -------
        object
            An object representing the trace of the marker in the scene.

        """
        marker_style = dict(size=size, color=color, symbol=marker_symbol)
        marker_trace = go.Scatter3d(
            x=position[0],
            y=position[1],
            z=position[2],
            marker=marker_style,
            name=label,
            showlegend=False if label is None else True,
        )
        self.figure.add_trace(marker_trace)
        return marker_trace

    def draw_sphere(self, position, *, color, label, radius):
        """Draws an sphere into the scene.

        Parameters
        ----------
        position : list[float, float]
            A list containing the x and y coordinates of the sphere location.
        color : str, optional
            A string representing the hexadecimal color for the sphere.
        label : str
            The name shown in the legend of the figure to identify the sphere.
        radius : float, optional
            The radius of the sphere.

        Returns
        -------
        object
            An object representing the trace of the sphere in the scene.

        """
        xx, yy, zz = generate_sphere(radius, position)
        sphere = go.Surface(
            x=xx,
            y=yy,
            z=zz,
            colorscale=[[0, color], [1, color]],
            cauto=False,
            cmin=1,
            cmax=1,
            showscale=False,
            name=label,
            showlegend=False if label is None else True,
        )
        self.figure.add_trace(sphere)
        return sphere

    def draw_coordinates(self, coordinates, *, colors, dashed, label):
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
            The name shown in the legend of the figure to identify the coordinates.

        Returns
        -------
        trace_coordinates : object
            An object representing the trace of the coordinates in the scene.

        """
        # Select the desired linestyle for the line representing the coordinates
        linestyle = "dash" if dashed else "solid"

        # Plot the coordinates in the scene
        coordinates_trace = go.Scatter3d(
            x=coordinates[0],
            y=coordinates[1],
            z=coordinates[2],
            line=dict(color=colors[0], width=5, dash=linestyle),
            mode="lines",
            name=label,
            showlegend=False if label is None else True,
        )
        self.figure.add_trace(coordinates_trace)
        return coordinates_trace

    def draw_axes_labels_with_length_scale_units(self, length_scale_units):
        """Draws the desired label into the specified axis.

        Parameters
        ----------
        lenght_scale_units : ~astropy.units.Unit
            Desired units of lenght used for representing distances.

        """
        self.figure.update_layout(
            scene=dict(
                xaxis=dict(title=f"x ({length_scale_units.name})"),
                yaxis=dict(title=f"y ({length_scale_units.name})"),
                zaxis=dict(title=f"z ({length_scale_units.name})"),
            )
        )
