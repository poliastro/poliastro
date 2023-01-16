"""A module implementing the base class for a orbit plotter backend."""


class OrbitPlotterBackend:
    """A base class for implementing new orbit plotter backends."""

    def __init__(self, scene, name):
        """Initialize the orbit plotter backend.

        Parameters
        ----------
        scene : object
            An instance representing the canvas or scene.
        name : str
            Name of the backend.

        Notes
        -----
        An orbit plotter backend instance gets initialized from a scene. This
        can be a :ref:`~matplotlib.Axes`, :ref:`~plotly.Figure` or any other
        object acting as canvas for rendering the scene.

        """
        # Verify backend name ends with '2D' or '3D'
        if name[-2:] not in ["2D", "3D"]:
            print(f"Name was found to be {name[-2:] = }")
            raise ValueError("Backend name must end with '2D' or '3D'.")

        self._scene = scene
        self._name = name

    @property
    def scene(self):
        """Return the scene object."""
        return self._scene

    @property
    def name(self):
        """Return the name of the backend.

        Returns
        -------
        str
            Name of the backend.

        """
        return self._name

    @property
    def is_2D(self):
        """Assert if backend is 2D.

        Returns
        -------
        bool
            ``True`` if it is a 2D backend, ``False`` if it is not.

        """
        return self.name.endswith("2D")

    @property
    def is_3D(self):
        """Assert if backend is 3D.

        Returns
        -------
        bool
            ``True`` if it is a 3D backend, ``False`` if it is not.

        """
        return self.name.endswith("3D")

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
        raise NotImplementedError(
            "This method is expected to be override by a plotting backend class."
        )

    def draw_marker(self, position, *, color, label, marker_symbol, size):
        """Draw desired marker into the scene.

        Parameters
        ----------
        position : list[float, float, float]
            A list containing the x, y and z coordinates of the point.
        color : str
            A string representing the hexadecimal color for the point.
        label : str
            The name to be used in the legend for the marker.
        marker_symbol : str
            The marker symbol to be used when drawing the point.
        size : float
            The size of the marker.

        """
        raise NotImplementedError(
            "This method is expected to be override by a plotting backend class."
        )

    def draw_position(self, position, *, color, label, size):
        """Draw the position of a body in the scene.

        Parameters
        ----------
        position : list[float, float, float]
            A list containing the x, y and z coordinates of the point.
        color : str
            A string representing the hexadecimal color for the marker.
        label : str
            The name to be used in the legend for the marker.
        size : float
            The size of the marker.

        Returns
        -------
        trace_position : object
            An object representing the trace of the position in the scene.

        """
        raise NotImplementedError(
            "This method is expected to be override by a plotting backend class."
        )

    def draw_impulse(self, position, *, color, label, size):
        """Draw an impulse into the scene.

        Parameters
        ----------
        position : list[float, float, float]
            A list containing the x, y and z coordinates of the impulse location.
        color : str
            A string representing the hexadecimal color for the impulse marker.
        label : str
            The name to be used in the legend for the marker.
        size : float
            The size of the marker for the impulse.

        """
        raise NotImplementedError(
            "This method is expected to be override by a plotting backend class."
        )

    def draw_sphere(self, position, *, color, label, radius):
        """Draw an sphere into the scene.

        Parameters
        ----------
        position : list[float, float, float]
            A list containing the x, y and z coordinates of the sphere location.
        color : str
            A string representing the hexadecimal color for the sphere.
        label : str
            The name to be used in the legend for the marker.
        radius : float
            The radius of the sphere.

        """
        raise NotImplementedError(
            "This method is expected to be override by a plotting backend class."
        )

    def undraw_attractor(self):
        """Remove the attractor from the scene."""
        raise NotImplementedError(
            "This method is expected to be override by a plotting backend class."
        )

    def draw_coordinates(self, coordinates, *, colors, label, size):
        """Draw desired coordinates into the scene.

        Parameters
        ----------
        position : list[list[float, float, float], ...]
            A set of lists containing the x, y and z coordinates of the sphere location.
        colors : list[str]
            A string representing the hexadecimal color for the coordinates.
        label : str
            The name to be used in the legend for the marker.
        size : float
            The size of the marker for drawing the coordinates.

        """
        raise NotImplementedError(
            "This method is expected to be override by a specific plotting backend."
        )

    def draw_axes_labels_with_length_scale_units(self, length_scale_units):
        """Draw the desired label into the specified axis.

        Parameters
        ----------
        lenght_scale_units : ~astropy.units.Unit
            Desired units of lenght used for representing distances.

        """
        raise NotImplementedError(
            "This method is expected to be override by a specific plotting backend."
        )

    def update_legend(self):
        """Update the legend of the scene."""
        raise NotImplementedError(
            "This method is expected to be override by a specific plotting backend."
        )

    def resize_limits(self):
        """Resize the limits of the scene."""
        raise NotImplementedError(
            "This method is expected to be override by a specific plotting backend."
        )

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

        """
        raise NotImplementedError(
            "This method is expected to be override by a specific plotting backend."
        )

    def show(self):
        """Display the scene."""
        raise NotImplementedError(
            "This method is expected to be override by a specific plotting backend."
        )
