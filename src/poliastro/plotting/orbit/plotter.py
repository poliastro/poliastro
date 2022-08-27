"""A module containing different orbit related plotters."""

from collections import namedtuple
import warnings

import astropy.units as u
import numpy as np

from poliastro.frames import Planes
from poliastro.plotting.orbit.backends import SUPPORTED_ORBIT_PLOTTER_BACKENDS
from poliastro.plotting.util import BODY_COLORS, generate_label
from poliastro.util import norm


class Trajectory(
    namedtuple(
        "Trajectory", ["coordinates", "position", "label", "colors", "dashed"]
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
        scene=None,
        backend_name="matplotlib2D",
        backend_kargs=None,
        num_points=150,
        *,
        plane=None,
        unit=u.km,
    ):
        """Initializes the plotter instance.

        Parameters
        ----------
        scene : object
            An instance representing the canvas or scene.
        backend_name : str
            Name of the plotting backend to be used.
        backend_kargs : dict, optional
            Additional configuration arguments for the backend.
        num_points : int, optional
            Number of points to use when drawing trajectories. Default to 150.
        plane : ~poliastro.frames.Plane, optional
            Reference plane to be used when drawing the scene. Default to
            `EARTH_EQUATOR`.
        unit : ~astropy.units.Unit
            Desired lenght unit to be used when representing distances.

        """
        # Verify selected backend is supported
        try:
            self._backend = SUPPORTED_ORBIT_PLOTTER_BACKENDS[backend_name](
                scene,
            )
        except KeyError:
            raise ValueError(
                f"Backend {desired_backend_name} is not supported. Available backends are {SUPPORTED_PLOTTING_BACKENDS.keys()}"
            )

        # Initialize the rest of the attributes
        self._attractor = None
        self._attractor_radius = np.inf * u.km
        self._num_points = num_points
        self._plane = plane or Planes.EARTH_EQUATOR
        self._frame = None
        self._trajectories = []  # type: List[Trajectory]
        self._unit = unit

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
        """Sets perifocal frame based on an orbit.

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
            elif not np.allclose([p_vec @ q_vec, q_vec @ w_vec, w_vec @ p_vec], 0):
                raise ValueError("Vectors must be mutually orthogonal.")
            else:
                self._frame = p_vec, q_vec, w_vec

            if self._trajectories:
                self._redraw()

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
        min_distance = min(
            [
                coordinates.norm().min()
                for coordinates, *_ in self._trajectories
            ]
            or [0 * u.m]
        )
        self._attractor_radius = max(
            self._attractor.R.to(u.km), min_distance.to(u.km) * 0.15
        )

        color = BODY_COLORS.get(self._attractor.name, "#999999")

        self._unplot_attractor()

        self.backend.draw_sphere(
            position=[0, 0, 0] * u.km,
            size=self._attractor_radius,
            color=color,
        )

    def _redraw(self):
        for trajectory in self._trajectories:
            self._plot_coordinates_and_position(trajectory)

    def _plot_position(self, position, label, colors):
        radius = min(
            self._attractor_radius * 0.5,
            (norm(position) - self._attractor.R) * 0.5,
        )  # Arbitrary thresholds
        self.backend.draw_point(radius, colors[0], label, center=position)

    def _plot_coordinates_and_position(self, trajectory):
        coordinates, position, label, colors, dashed = trajectory

        trace_coordinates = self._plot_coordinates(
            coordinates, label, colors, dashed
        )

        if position is not None:
            trace_position = self._plot_position(position, label, colors)
        else:
            trace_position = None

        return trace_coordinates, trace_position

    def _add_trajectory(
        self, coordinates, position=None, *, colors=None, dashed=False, label=None,
    ):
        """Add a new trajectory to the scene.

        Parameters
        ----------
        coordinates : ~astropy.coordinates.CartesianRepresentation
            Trajectory to plot.
        position : ~numpy.ndarray, optional
            Position of the body along its orbit.
        colors : str, optional
            A list of colors for the points of the trajectory.
        dashed : bool, optional
            Fade the orbit trail, default to False.
        label : str, optional
            Label of the orbit.

        Returns
        -------
        trace_coordinates : object
            An object representing the trace of the coordinates in the scene.
        trace_position : object
            An object representing the trace of the position in the scene.

        """
        # Instantiate the trajectory and append it to the internal list
        trajectory = Trajectory(coordinates, position, label, colors, dashed)
        self._trajectories.append(trajectory)

        # Update the attractor appearance (if required) based on new scene trajectory
        self._plot_attractor()

        # Project the coordinates into desired frame for 2D backends
        if self.backend.is_2D:
            rr = coordinates.xyz.transpose()
            coordinates = self._project(rr)

        # Add the coordinates and the position to the scene
        trace_coordinates = self.backend.draw_coordinates(coordinates, color=colors)
        trace_position = self.backend.draw_position(position, color=colors)

        return trace_coordinates, trace_position

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
        # Check if the orbit plotter frame has been set
        if not self._frame:
            self.set_orbit_frame(orbit)

        # Assign attractor if required
        self.set_attractor(orbit.attractor)

        # Represent the orbit w.r.t. the derired plane
        orbit = orbit.change_plane(self.plane)

        # Get orbit colors and label
        colors = self.backend._get_colors(color, trail)
        label = generate_label(orbit.epoch, label)

        # Compute the coorindates representing the trajectory
        coordinates = orbit.sample(self._num_points)

        return self._add_trajectory(
            coordinates, orbit.r, label=label, colors=colors, dashed=True
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
        self.set_attractor(body.parent)

        if color is None:
            color = BODY_COLORS.get(body.name)

        # Get approximate, mean value for the period
        period = get_mean_elements(body, epoch).period

        label = generate_label(epoch, label or str(body))
        epochs = time_range(
            epoch, num_values=self._num_points, end=epoch + period, scale="tdb"
        )
        ephem = Ephem.from_body(
            body, epochs, attractor=body.parent, plane=self.plane
        )

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
        r0 = ephem.rv(epoch)[0] if epoch else None

        # Get the list of colors for previous coordinate points
        colors = self.backend._get_colors(color, trail)

        return self._add_trajectory(
            coordinates, r0, label=str(label), colors=colors, dashed=False
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

        if len(maneuver_phases) == 0:
            # For single-impulse maneuver only draw the impulse marker
            impulse_label = f"Impulse 1 - {label}"
            impulse_lines = (
                [self._draw_impulse(color, impulse_label, final_phase.r)],
            )
            return [(impulse_label, impulse_lines)]
        else:
            # Declare for holding (label, lines) for each impulse and trajectory
            lines_list = []

            # Collect the coordinates for the different maneuver phases
            for ith_impulse, orbit_phase in enumerate(maneuver_phases):

                # Plot the impulse marker and collect its label and lines
                impulse_label = f"Impulse {ith_impulse + 1} - {label}"
                impulse_lines = (
                    [self._draw_impulse(color, impulse_label, orbit_phase.r)],
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
                trajectory_lines = self._plot_trajectory(
                    phase_coordinates, label=label, color=color, trail=trail
                )
                lines_list.append((label, trajectory_lines))

            # Finally, draw the impulse at the very beginning of the final phase
            impulse_label = f"Impulse {ith_impulse + 2} - {label}"
            impulse_lines = (
                [self._draw_impulse(color, impulse_label, final_phase.r)],
            )
            lines_list.append((impulse_label, impulse_lines))

    def plot_trajectory(
        self, coordinates, *, label=None, color=None, trail=False
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

        Raises
        ------
        ValueError
            An attractor must be set first.
        """
        if self._attractor is None:
            raise ValueError(
                "An attractor must be set up first, please use "
                "set_attractor(Major_Body) or plot(orbit)"
            )

        colors = self._get_colors(color, trail)
        coordinates = coordinates.represent_as(CartesianRepresentation)

        return self._add_trajectory(
            coordinates, None, label=str(label), colors=colors, dashed=False
        )

    def show(self):
        """Render the plot."""
        self._backend.show()
