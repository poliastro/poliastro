"""A module containing different orbit related plotters."""

import warnings
from collections import namedtuple
from typing import List

import astropy.units as u
import numpy as np
from astropy.coordinates import CartesianRepresentation

import poliastro.plotting.orbit.backends as orbit_plotter_backends
from poliastro.ephem import Ephem
from poliastro.frames import Planes
from poliastro.plotting.util import BODY_COLORS, generate_label
from poliastro.twobody.mean_elements import get_mean_elements
from poliastro.twobody.sampling import EpochBounds
from poliastro.util import norm, time_range


class Trajectory(
    namedtuple(
        "Trajectory", ["coordinates", "position", "colors", "dashed", "label"]
    )
):
    """A class for collecting all information of a body within a plotter.

    Information contains all the trajectory coordinates, the current position,
    the name of the label to be displayed in the plotter, the color to be used
    and if dashed mode is desired or not.

    """

    pass


class OrbitPlotter:
    """A base class containing common attributes and methods for plotters."""

    def __init__(
        self,
        backend=None,
        num_points=150,
        *,
        plane=None,
        length_scale_units=u.km,
    ):
        """Initializes the plotter instance.

        Parameters
        ----------
        backend : ~poliastro.plotting.orbit.backends._base.OrbitPlotterBackend
            An instance of ``OrbitPlotterBackend`` for rendendering the scene.
        num_points : int, optional
            Number of points to use when drawing trajectories. Default to 150.
        plane : ~poliastro.frames.Plane, optional
            Reference plane to be used when drawing the scene. Default to
            `EARTH_EQUATOR`.
        length_scale_units : ~astropy.units.Unit
            Desired units of length used for representing distances.
        _
        """
        # Initialize the backend, number of points and length scale
        self._backend = backend or orbit_plotter_backends.Matplotlib2D(
            ax=None, use_dark_theme=False
        )
        self._num_points = num_points
        self._length_scale_units = length_scale_units

        # Initialize the attractor
        self._attractor = None
        self._attractor_radius = np.inf * length_scale_units

        # Initialize the plane and frame used by the plotter
        self._plane = plane or Planes.EARTH_EQUATOR
        self._frame = None

        # Initialize the list containing all the plotted trajectories
        self._trajectories = []  # type: List[Trajectory]

    @property
    def backend(self):
        """Backend instance used by the plotter."""
        return self._backend

    @property
    def plane(self):
        """Reference plane to be used when drawing the scene."""
        return self._plane

    @property
    def trajectories(self):
        """A list with all the `Trajectory` instances used in the plotter."""
        return self._trajectories

    @property
    def length_scale_units(self):
        """Return the units of length used for representing distances.

        Returns
        -------
        length_units : ~astropy.units.Unit
            Desired length units to be used when representing distances.

        """
        return self._length_scale_units

    def set_attractor(self, attractor):
        """Set the desired plotting attractor.

        Parameters
        ----------
        attractor : ~poliastro.bodies.Body
            Central body.

        Raises
        ------
        NotImplementedError
            Raised if attractor is already set.

        """
        if self._attractor is None:
            self._attractor = attractor
        elif attractor is not self._attractor:
            raise NotImplementedError(
                f"Attractor has already been set to {self._attractor.name}"
            )

    def set_orbit_frame(self, orbit):
        """Set the perifocal frame based on an orbit.

        Parameters
        ----------
        orbit : ~poliastro.twobody.Orbit
            Orbit to use as frame.

        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)

            p_vec, q_vec, w_vec = orbit.pqw()

            if not np.allclose([norm(v) for v in (p_vec, q_vec, w_vec)], 1):
                raise ValueError("Vectors must be unit.")
            elif not np.allclose(
                [p_vec @ q_vec, q_vec @ w_vec, w_vec @ p_vec], 0
            ):
                raise ValueError("Vectors must be mutually orthogonal.")
            else:
                self._frame = p_vec, q_vec, w_vec

            if self._trajectories:
                self._redraw()

    def set_body_frame(self, body, epoch=None):
        """Sets perifocal frame based on the orbit of a body at a particular epoch if given.

        Parameters
        ----------
        body : poliastro.bodies.SolarSystemPlanet
            Body.
        epoch : astropy.time.Time, optional
            Epoch of current position.

        """
        from warnings import warn

        from astropy import time

        # HACK: avoid circular dependency with ``Body.plot()``
        from poliastro.bodies import Sun
        from poliastro.twobody import Orbit
        from poliastro.warnings import TimeScaleWarning

        if not epoch:
            epoch = time.Time.now().tdb
        elif epoch.scale != "tdb":
            epoch = epoch.tdb
            warn(
                "Input time was converted to scale='tdb' with value "
                f"{epoch.tdb.value}. Use Time(..., scale='tdb') instead.",
                TimeScaleWarning,
                stacklevel=2,
            )

        with warnings.catch_warnings():
            ephem = Ephem.from_body(body, epoch, attractor=Sun, plane=self.plane)  # type: ignore
            orbit = Orbit.from_ephem(Sun, ephem, epoch).change_plane(self.plane)  # type: ignore

        self.set_orbit_frame(orbit)

    def _project(self, vec):
        """Project the vector into the frame of the orbit plotter.

        Parameters
        ----------
        vec : ~np.ndarray
            The vector to be projected into the frame of the orbit plotter.

        Returns
        -------
        tuple(float, float, float)
            A tuple containing the x, y, and z coordinates of the vector
            projected into the frame of the orbit plotter.

        """
        vec_proj = vec - (vec @ self._frame[2])[:, None] * self._frame[2]
        return [vec_proj @ self._frame[i] for i in range(3)]

    def _unplot_attractor(self):
        """Removes the attractor from the scene."""
        self._backend.undraw_attractor()

    def _plot_attractor(self):
        """Plot the scene attractor.

        Notes
        -----
        This method allows to select a sensible value for the radius of the
        attractor.

        """
        ref_len_units = self.length_scale_units
        min_distance = min(
            [
                coordinates.norm().min().to(ref_len_units)
                for coordinates, *_ in self._trajectories
            ]
            or [0 * ref_len_units]
        )
        self._attractor_radius = max(
            self._attractor.R.to(ref_len_units),
            min_distance.to(ref_len_units) * 0.15,
        )

        color = BODY_COLORS.get(self._attractor.name, "#999999")

        self._unplot_attractor()

        self.backend.draw_sphere(
            position=[0, 0, 0],
            color=color,
            label=None,
            radius=self._attractor_radius.to_value(ref_len_units),
        )

    def _redraw(self):
        for trajectory in self._trajectories:
            self._add_trajectory(trajectory)

    def _create_trajectory(
        self, coordinates, position, *, colors=None, dashed=False, label=None
    ):
        """Create a new ``Trajectory`` instance.

        Parameters
        ----------
        coordinates : ~astropy.coordinates.CartesianRepresentation
            Trajectory to plot.
        position : ~numpy.ndarray, optional
            Position of the body along its orbit.
        colors : str, optional
            A list of colors for the points of the trajectory.
        dashed : bool, optional
            ``True`` to use a dashed line style. ``False`` otherwise.
        label : str, optional
            Label of the orbit.

        Returns
        -------
        Trajectory
            An object for modeling trajectories.

        """
        # Instantiate the trajectory and append it to the internal list
        trajectory = Trajectory(coordinates, position, colors, dashed, label)
        self._trajectories.append(trajectory)
        return trajectory

    def _add_trajectory(self, trajectory):
        """Add a new trajectory to the scene.

        Parameters
        ----------
        trajectory : Trajectory
            An object for modeling trajectories.

        Returns
        -------
        trace_coordinates : object
            An object representing the trace of the coordinates in the scene.
        trace_position : object
            An object representing the trace of the position in the scene.

        """
        # Update the attractor appearance (if required) based on new scene trajectory
        self._plot_attractor()

        # Unpack trajectory data
        coordinates, position, colors, dashed, label = trajectory

        # HACK: the different backends manage legends in a different way. There
        # are additional restrictions depending on the type of objects they
        # provide for generating the coordinates and the positions. Some of
        # these objects can not be rendered in the legend.
        has_coordinates, has_position = (
            coordinates is not None,
            position is not None,
        )
        coordinates_label, position_label = self.backend.generate_labels(
            label, has_coordinates, has_position
        )

        # Project the coordinates into desired frame for 2D backends
        if self.backend.is_2D:
            rr = coordinates.xyz.transpose()
            coordinates = self._project(rr)
            if position is not None:
                position = (
                    np.asarray(
                        self._project(
                            [position.to_value(self.length_scale_units)]
                        )
                    ).flatten()
                    * self.length_scale_units
                )
        else:
            coordinates = (
                coordinates.x,
                coordinates.y,
                coordinates.z,
            )

        # Add the coordinates to the scene
        trace_coordinates = self.backend.draw_coordinates(
            [
                coords.to_value(self.length_scale_units)
                for coords in coordinates
            ],
            colors=colors,
            dashed=dashed,
            label=coordinates_label,
        )

        # Add the position to the scene
        if position is not None:
            # Use a proper size for the position in the scene
            marker_size = min(
                self._attractor_radius * 0.5,
                (norm(position) - self._attractor.R) * 0.5,
            )

            # Generate the trace for the position
            trace_position = self.backend.draw_position(
                position.to_value(self.length_scale_units),
                color=colors[0],
                label=position_label,
                size=marker_size.to_value(self.length_scale_units),
            )
        else:
            trace_position = None

        # Update the legend to reflect the new traces
        if coordinates_label is not None or position_label is not None:
            self.backend.update_legend()

        # Update the limits of the scene
        self.backend.resize_limits()

        # Update the axis legends using the desired length scale units
        self.backend.draw_axes_labels_with_length_scale_units(
            self.length_scale_units
        )

        return (
            (trace_coordinates, trace_position)
            if position is not None
            else (trace_coordinates, None)
        )

    def plot(self, orbit, *, color=None, label=None, trail=False, dashed=True):
        """Plots state and osculating orbit in their plane.

        Parameters
        ----------
        orbit : ~poliastro.twobody.orbit.Orbit
            Orbit to plot.
        color : str, optional
            Color of the line and the position.
        label : str, optional
            Label of the orbit.
        trail : bool, optional
            Fade the orbit trail, default to False.
        dashed : bool, optional
            ``True`` to use a dashed line style. ``False`` otherwise.

        """
        # Check if the orbit plotter frame has been set
        if self.backend.is_2D and self._frame is None:
            self.set_orbit_frame(orbit)

        # Assign attractor if required
        self.set_attractor(orbit.attractor)

        # Represent the orbit w.r.t. the derired plane
        orbit = orbit.change_plane(self.plane)

        # Get orbit label
        label = generate_label(orbit.epoch, label)

        # Compute the coordinates and body position alongs its orbit
        coordinates = orbit.sample(self._num_points)
        position = orbit.r

        return self.plot_coordinates(
            coordinates,
            position=position,
            label=label,
            color=color,
            trail=trail,
            dashed=dashed,
        )

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
        # Assign frame if required
        if self.backend.is_2D and self._frame is None:
            self.set_body_frame(body, epoch)

        # Assign attractor if required
        self.set_attractor(body.parent)

        # Get approximate, mean value for the period and generate ephemerides
        period = get_mean_elements(body, epoch).period
        epochs = time_range(
            epoch, num_values=self._num_points, end=epoch + period, scale="tdb"
        )
        ephem = Ephem.from_body(
            body, epochs, attractor=body.parent, plane=self.plane
        )

        # Get body color and label
        if color is None:
            color = BODY_COLORS.get(body.name)
        label = generate_label(epoch, label or str(body))

        return self.plot_ephem(
            ephem, epoch, label=label, color=color, trail=trail
        )

    def plot_ephem(
        self, ephem, epoch=None, *, label=None, color=None, trail=False
    ):
        """Plots Ephem object over its sampling period.

        Parameters
        ----------
        ephem : ~poliastro.ephem.Ephem
            Ephemerides to plot.
        epoch : astropy.time.Time, optional
            Epoch of the current position, `None` is used if not given.
        label : str, optional
            Label of the orbit, default to the name of the body.
        color : str, optional
            Color of the line and the position.
        trail : bool, optional
            Fade the orbit trail, default to False.

        """
        if self.backend.is_2D and self._frame is None:
            raise ValueError(
                "A frame must be set up first, please use "
                "set_orbit_frame(orbit) or plot(orbit)"
            )

        if self._attractor is None:
            raise ValueError(
                "An attractor must be set up first, please use "
                "set_attractor(MajorBody) or plot(orbit)"
            )

        if ephem.plane is not self.plane:
            raise ValueError(
                f"The ephemerides reference plane is {ephem.plane} "
                f"while the plotter is using {self.plane}, "
                "sample the ephemerides using a different plane "
                "or create a new plotter"
            )

        # Collect the coordinates of the trajectory defined by the ephemerides
        coordinates = ephem.sample()

        # Get the location in the epehemerides associated with the desired epoch
        position = ephem.rv(epoch)[0] if epoch is not None else None

        return self.plot_coordinates(
            coordinates,
            position=position,
            label=str(label),
            color=color,
            trail=trail,
        )

    def plot_maneuver(
        self, initial_orbit, maneuver, label=None, color=None, trail=False
    ):
        """Plots the maneuver trajectory applied to the provided initial orbit.

        Parameters
        ----------
        initial_orbit : ~poliastro.twobody.orbit.Orbit
            The base orbit for which the maneuver will be applied.
        maneuver : ~poliastro.maneuver.Maneuver
            The maneuver to be plotted.
        label : str, optional
            Label of the trajectory.
        color : str, optional
            Color of the trajectory.
        trail : bool, optional
            Fade the orbit trail, default to False.

        """
        if self._attractor is None:
            raise ValueError(
                "An attractor must be set up first, please use "
                "set_attractor(Major_Body) or plot(orbit)"
            )

        # Apply the maneuver, collect all intermediate states and allocate the
        # final coordinates list array
        *maneuver_phases, final_phase = initial_orbit.apply_maneuver(
            maneuver, intermediate=True
        )

        # Project the coordinates into desired frame for 2D backends
        if self.backend.is_2D:
            final_phase_position = (
                np.asarray(
                    self._project(
                        [final_phase.r.to_value(self.length_scale_units)]
                    )
                ).flatten()
                * self.length_scale_units
            )
        else:
            final_phase_position = final_phase.r

        if not len(maneuver_phases):
            # For single-impulse maneuver only draw the impulse marker
            impulse_label = f"Impulse - {label}"
            impulse_lines = (
                [
                    self.backend.draw_impulse(
                        position=final_phase_position.to_value(
                            self.length_scale_units
                        ),
                        color=color,
                        label=impulse_label,
                        size=None,
                    )
                ],
            )
            return [(impulse_label, impulse_lines)]
        else:
            # Declare for holding (label, lines) for each impulse and trajectory
            lines_list = []

            # Collect the coordinates for the different maneuver phases
            for ith_impulse, orbit_phase in enumerate(maneuver_phases):

                # Project the coordinates into desired frame for 2D backends
                if self.backend.is_2D:
                    orbit_phase_r = (
                        np.asarray(
                            self._project(
                                [
                                    orbit_phase.r.to_value(
                                        self.length_scale_units
                                    )
                                ]
                            )
                        ).flatten()
                        * self.length_scale_units
                    )
                else:
                    orbit_phase_r = orbit_phase.r

                # Plot the impulse marker and collect its label and lines
                impulse_label = f"Impulse {ith_impulse + 1} - {label}"
                impulse_lines = self.backend.draw_impulse(
                    position=orbit_phase_r.to_value(self.length_scale_units),
                    color=color,
                    label=impulse_label,
                    size=None,
                )
                lines_list.append((impulse_label, impulse_lines))

                # HACK: if no color is provided, get the one randomly generated
                # for previous impulse lines
                color = (
                    impulse_lines[0][0].get_color() if color is None else color
                )

                # Get the propagation time required before next impulse
                time_to_next_impulse, _ = maneuver.impulses[ith_impulse + 1]

                # Collect the coordinate points for the i-th orbit phase
                # TODO: Remove `.sample()` to return Ephem and use `plot_ephem` instead?
                phase_coordinates = orbit_phase.to_ephem(
                    strategy=EpochBounds(
                        min_epoch=orbit_phase.epoch,
                        max_epoch=orbit_phase.epoch + time_to_next_impulse,
                    )
                ).sample()

                # Plot the phase trajectory and collect its label and lines
                trajectory_lines = self.plot_coordinates(
                    phase_coordinates,
                    position=None,
                    color=color,
                    label=label,
                    trail=trail,
                )
                lines_list.append((label, trajectory_lines))

            # Finally, draw the impulse at the very beginning of the final phase
            impulse_label = f"Impulse {ith_impulse + 2} - {label}"
            impulse_lines = self.backend.draw_impulse(
                position=final_phase_position.to_value(
                    self.length_scale_units
                ),
                color=color,
                label=impulse_label,
                size=None,
            )
            lines_list.append((impulse_label, impulse_lines))

            # Update the legend to reflect all the traces
            self.backend.update_legend()

    def plot_coordinates(
        self,
        coordinates,
        *,
        position=None,
        label=None,
        color=None,
        trail=False,
        dashed=False,
    ):
        """Plots a precomputed trajectory.

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
        dashed : bool, optional
            ``True`` to use a dashed line style. ``False`` otherwise.

        Raises
        ------
        ValueError
            An attractor must be set first.

        """
        # Check if a frame has been set
        if self.backend.is_2D and self._frame is None:
            raise ValueError(
                "A frame must be set up first, please use "
                "set_orbit_frame(orbit) or plot(orbit)"
            )

        # Check if the attractor
        if self._attractor is None:
            raise ValueError(
                "An attractor must be set up first, please use "
                "set_attractor(Major_Body) or plot(orbit)"
            )

        # Get orbit colors and label
        colors = self.backend._get_colors(color=color, trail=trail)

        # Force Cartesian representation for coordinates
        coordinates = coordinates.represent_as(CartesianRepresentation)

        # Generate the trajectory instance
        trajectory = self._create_trajectory(
            coordinates,
            position,
            colors=colors,
            dashed=dashed,
            label=label,
        )
        return self._add_trajectory(trajectory)

    def plot_trajectory(
        self,
        coordinates,
        *,
        label=None,
        color=None,
        trail=False,
        dashed=False,
    ):
        """Plots a precomputed trajectory.

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
        dashed : bool, optional
            ``True`` to use a dashed line style. ``False`` otherwise.

        Raises
        ------
        ValueError
            An attractor must be set first.

        """
        return self.plot_coordinates(
            coordinates,
            position=None,
            color=color,
            label=label,
            trail=trail,
            dashed=dashed,
        )

    @u.quantity_input(elev=u.rad, azim=u.rad, distance=u.km)
    def set_view(self, elevation_angle, azimuth_angle, distance=5 * u.km):
        """Changes 3D view by setting the elevation, azimuth and distance.

        Parameters
        ----------
        elevation_angle : ~astropy.units.Quantity
            Desired elevation angle of the camera.
        azimuth_angle : ~astropy.units.Quantity
            Desired azimuth angle of the camera.
        distance : optional, ~astropy.units.Quantity
            Desired distance of the camera to the scene.

        """
        if self.backend.is_2D:
            raise AttributeError("View can only be in 3D backends.")
        self.backend.set_view(
            elevation_angle.to_value(u.rad),
            azimuth_angle.to_value(u.rad),
            distance.to_value(self.length_scale_units),
        )

    def show(self):
        """Render the plot."""
        self.backend.show()
