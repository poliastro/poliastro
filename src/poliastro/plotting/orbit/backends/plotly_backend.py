"""A module implementing orbit plotter backends based on Plotly."""

from itertools import cycle

import numpy as np
import plotly
import plotly.graph_objects as go
from astropy import units as u

from poliastro.plotting.orbit.backends.base_backend import _OrbitPlotterBackend
from poliastro.plotting.util import generate_sphere


class OrbitPlotterBackendPlotly(_OrbitPlotterBackend):
    """An orbit plotter backend class based on Plotly."""

    def __init__(self, figure=None, layout=None):
        """Initializes a backend instance.

        Parameters
        ----------
        figure : ~plotly.graph_objects.Figure
            The plotly ``Figure`` to render the scene.
        ~plotly.graph_objects.Layout
            The plotly ``Layout`` object linked to the figure.

        """
        figure = figure or go.Figure()
        super().__init__(figure, name=self.__class__.__name__)

        self._layout = layout or go.Layout()
        self.update_layout()

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

    def update_layout(self, layout=None):
        """Update the layout of the figure scene.

        Parameters
        ----------
        layout : ~plotly.graph_objects.Layout
            The new plotly ``Layout`` to be used in the figure.

        """
        layout = self.layout or layout
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
        # This will apply the label to either the point or the osculating
        # orbit depending on the last plotted line
        trace_id = -1 if trace_position is None else -2
        self.figure["data"][trace_id]["showlegend"] = True
        self.figure["data"][trace_id]["name"] = label

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
        return self.draw_sphere(position, color=color, radius=size)

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

    def show(self):
        """Displays the scene."""
        self.update_layout()
        if not self.figure._in_batch_mode:
            return self.figure.show()


class OrbitPlotterBackendPlotly2D(OrbitPlotterBackendPlotly):
    """An orbit plotter backend class based on Plotly."""

    def __init__(self, figure=None, use_dark_theme=None):
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
            xaxis=dict(title=" x ({self._unit})", constrain="domain"),
            yaxis=dict(title=" y ({self._unit})", scaleanchor="x"),
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
            A string containing tha message for the label.
        size : float, optional
            Desired size for the marker.

        Returns
        -------
        ~matplotlib.lines.Line2D
            An object representing the trace of the marker in the scene.

        """
        marker = dict(size=size, color=color, symbol=marker_symbol)
        marker_trace = go.Scatter(
            x=position[0],
            y=position[1],
            marker=marker,
            name=label,
            showlegend=True if label is not None else False,
        )
        self.figure.add_trace(marker_trace)
        return marker_trace

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
        return self.figure.add_shape(
            type="circle",
            xref="x",
            yref="y",
            x0=(position[0] - radius).to_value(u.km),
            y0=(position[1] - radius).to_value(u.km),
            x1=(position[0] + radius).to_value(u.km),
            y1=(position[1] + radius).to_value(u.km),
            fillcolor=color,
            line_color=color,
            opacity=1,
        )

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

        Returns
        -------
        trace_coordinates : ~matplotlib.lines.Line2D
            An object representing the trace of the coordinates in the scene.

        """
        # Select the desired linestyle for the line representing the coordinates
        linestyle = "dash" if dashed else "solid"

        # Unpack coordinates
        x, y, _ = (coords.to_value(u.km) for coords in coordinates)

        # Plot the coordinates in the scene
        coordinates_trace = go.Scatter(
            x=x,
            y=y,
            line=dict(color=colors[0], width=5, dash=linestyle),
            mode="lines",
            showlegend=False,
        )
        self.figure.add_trace(coordinates_trace)
        return coordinates_trace


class OrbitPlotterBackendPlotly3D(OrbitPlotterBackendPlotly):
    """An orbit plotter backend class based on Plotly."""

    def __init__(self, figure=None, use_dark_theme=None):
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
                xaxis=dict(title=" x ({self._unit})"),
                yaxis=dict(title=" y ({self._unit})"),
                zaxis=dict(title=" z ({self._unit})"),
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
            A string containing tha message for the label.
        size : float, optional
            Desired size for the marker.

        Returns
        -------
        ~matplotlib.lines.Line2D
            An object representing the trace of the marker in the scene.

        """
        marker = dict(size=size, color=color, symbol=marker_symbol)
        marker_trace = go.Scatter3d(
            x=position[0],
            y=position[1],
            z=position[2],
            marker=marker,
            name=label,
            showlegend=True if label is not None else False,
        )
        self.figure.add_trace(marker_trace)
        return marker_trace

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
        xx, yy, zz = generate_sphere(radius, position)
        sphere = go.Surface(
            x=xx.to_value(u.km),
            y=yy.to_value(u.km),
            z=zz.to_value(u.km),
            colorscale=[[0, color], [1, color]],
            cauto=False,
            cmin=1,
            cmax=1,
            showscale=False,
        )
        self.figure.add_trace(sphere)
        return sphere

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

        Returns
        -------
        trace_coordinates : ~matplotlib.lines.Line2D
            An object representing the trace of the coordinates in the scene.

        """
        # Select the desired linestyle for the line representing the coordinates
        linestyle = "dash" if dashed else "solid"

        # Plot the coordinates in the scene
        coordinates_trace = go.Scatter3d(
            x=coordinates.x.to_value(u.km),
            y=coordinates.y.to_value(u.km),
            z=coordinates.z.to_value(u.km),
            line=dict(color=colors[0], width=5, dash=linestyle),
            mode="lines",
            showlegend=False,
        )
        self.figure.add_trace(coordinates_trace)
        return coordinates_trace

    @u.quantity_input(elev=u.rad, azim=u.rad, distance=u.km)
    def set_view(self, elev, azim, distance=5 * u.km):
        """Changes 3D view by setting the elevation, azimuth and distance.

        Parameters
        ----------
        elev : ~astropy.units.Quantity
            Desired elevation angle of the camera.
        azim : ~astropy.units.Quantity
            Desired azimuth angle of the camera.
        distance : optional, ~astropy.units.Quantity
            Desired distance of the camera to the scene.

        """
        x = distance * np.cos(elev) * np.cos(azim)
        y = distance * np.cos(elev) * np.sin(azim)
        z = distance * np.sin(elev)

        self.layout.update(
            {
                "scene": {
                    "camera": {
                        "eye": {
                            "x": x.to_value(u.km),
                            "y": y.to_value(u.km),
                            "z": z.to_value(u.km),
                        }
                    }
                }
            }
        )
