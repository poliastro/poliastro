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

    def __init__(self, figure, layout, ref_units):
        """Initializes a backend instance.

        Parameters
        ----------
        figure : ~plotly.graph_objects.Figure
            The plotly ``Figure`` to render the scene.
        ~plotly.graph_objects.Layout
            The plotly ``Layout`` object linked to the figure.
        ref_units : ~astropy.units.Unit
            Desired lenght units to be used when representing distances.

        """
        figure = figure or go.Figure()
        super().__init__(figure, self.__class__.__name__, ref_units)

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

    def show(self):
        """Displays the scene."""
        self.update_layout(self._layout)
        if not self.figure._in_batch_mode:
            return self.figure.show()


class OrbitPlotterBackendPlotly2D(OrbitPlotterBackendPlotly):
    """An orbit plotter backend class based on Plotly."""

    def __init__(self, figure, use_dark_theme, ref_units):
        """Initializes a backend instance.

        Parameters
        ----------
        figure : ~plotly.graph_objects.Figure
            The plotly ``Figure`` to render the scene.
        use_dark_theme : bool, optional
            If ``True``, uses dark theme. If ``False``, uses light theme.
            Default to ``False``.
        ref_units : ~astropy.units.Unit
            Desired lenght units to be used when representing distances.

        """
        # Apply the desired theme
        theme = "plotly_dark" if use_dark_theme is True else "plotly"

        # Declare the layout and attach it to the figure
        layout = go.Layout(
            autosize=True,
            xaxis=dict(title=f" x ({ref_units.name})", constrain="domain"),
            yaxis=dict(title=f" y ({ref_units.name})", scaleanchor="x"),
            template=theme,
        )
        super().__init__(figure, layout, ref_units)

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
            x0=(position[0] - radius).to_value(self.ref_units),
            y0=(position[1] - radius).to_value(self.ref_units),
            x1=(position[0] + radius).to_value(self.ref_units),
            y1=(position[1] + radius).to_value(self.ref_units),
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
        x, y, _ = (coords.to_value(self.ref_units) for coords in coordinates)

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

    def generate_labels(self, label, has_coordinates, has_position):
        return (label, None)


class OrbitPlotterBackendPlotly3D(OrbitPlotterBackendPlotly):
    """An orbit plotter backend class based on Plotly."""

    def __init__(self, figure, use_dark_theme, ref_units):
        """Initializes a backend instance.

        Parameters
        ----------
        figure : ~plotly.graph_objects.Figure
            The plotly ``Figure`` to render the scene.
        use_dark_theme : bool, optional
            If ``True``, uses dark theme. If ``False``, uses light theme.
            Default to ``False``.
        ref_units : ~astropy.units.Unit
            Desired lenght unit to be used when representing distances.

        """
        # Apply the desired theme
        theme = "plotly_dark" if use_dark_theme is True else "plotly"

        # Declare the layout and attach it to the figure
        layout = go.Layout(
            autosize=True,
            scene=dict(
                xaxis=dict(title=f" x ({ref_units.name})"),
                yaxis=dict(title=f" y ({ref_units.name})"),
                zaxis=dict(title=f" z ({ref_units.name})"),
                aspectmode="data",
            ),
            template=theme,
        )
        super().__init__(figure, layout, ref_units)

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
            x=xx.to_value(self.ref_units),
            y=yy.to_value(self.ref_units),
            z=zz.to_value(self.ref_units),
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
            x=coordinates.x.to_value(self.ref_units),
            y=coordinates.y.to_value(self.ref_units),
            z=coordinates.z.to_value(self.ref_units),
            line=dict(color=colors[0], width=5, dash=linestyle),
            mode="lines",
            name=label,
            showlegend=False if label is None else True,
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
                            "x": x.to_value(self.ref_units),
                            "y": y.to_value(self.ref_units),
                            "z": z.to_value(u.km),
                        }
                    }
                }
            }
        )

    def generate_labels(self, label, has_coordinates, has_position):
        return (label, None)
