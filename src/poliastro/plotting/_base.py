import warnings
from collections import namedtuple
from typing import List

import numpy as np
from astropy import units as u
from astropy.coordinates import CartesianRepresentation

from ..ephem import Ephem
from ..frames import Planes
from ..twobody.mean_elements import get_mean_elements
from ..util import norm, time_range
from .util import BODY_COLORS, generate_label


class Trajectory(
    namedtuple("Trajectory", ["coordinates", "position", "label", "colors", "dashed"])
):
    pass


class BaseOrbitPlotter:
    """
    Base class for all the OrbitPlotter classes.
    """

    def __init__(self, num_points=150, *, plane=None):
        self._num_points = num_points
        self._trajectories = []  # type: List[Trajectory]
        self._attractor = None
        self._plane = plane or Planes.EARTH_EQUATOR
        self._attractor_radius = np.inf * u.km

    @property
    def trajectories(self):
        return self._trajectories

    @property
    def plane(self):
        return self._plane

    def _set_attractor(self, attractor):
        if self._attractor is None:
            self._attractor = attractor
        elif attractor is not self._attractor:
            raise NotImplementedError(
                f"Attractor has already been set to {self._attractor.name}"
            )

    def set_attractor(self, attractor):
        """Sets plotting attractor.

        Parameters
        ----------
        attractor : ~poliastro.bodies.Body
            Central body.

        """
        self._set_attractor(attractor)

    def _clear_attractor(self):
        raise NotImplementedError

    def _redraw_attractor(self):
        # Select a sensible value for the radius: realistic for low orbits,
        # visible for high and very high orbits
        min_distance = min(
            [coordinates.norm().min() for coordinates, *_ in self._trajectories]
            or [0 * u.m]
        )
        self._attractor_radius = max(
            self._attractor.R.to(u.km), min_distance.to(u.km) * 0.15
        )

        color = BODY_COLORS.get(self._attractor.name, "#999999")

        self._clear_attractor()

        self._draw_sphere(
            self._attractor_radius,
            color,
            self._attractor.name,
        )

    def _redraw(self):
        for trajectory in self._trajectories:
            self.__plot_coordinates_and_position(trajectory)

    def _get_colors(self, color, trail):
        raise NotImplementedError

    def _draw_point(self, radius, color, name, center=None):
        raise NotImplementedError

    def _draw_sphere(self, radius, color, name, center=None):
        raise NotImplementedError

    def _plot_coordinates(self, coordinates, label, colors, dashed):
        raise NotImplementedError

    def _plot_position(self, position, label, colors):
        radius = min(
            self._attractor_radius * 0.5, (norm(position) - self._attractor.R) * 0.5
        )  # Arbitrary thresholds
        self._draw_point(radius, colors[0], label, center=position)

    def __plot_coordinates_and_position(self, trajectory):
        coordinates, position, label, colors, dashed = trajectory

        trace_coordinates = self._plot_coordinates(coordinates, label, colors, dashed)

        if position is not None:
            trace_position = self._plot_position(position, label, colors)
        else:
            trace_position = None

        return trace_coordinates, trace_position

    def __add_trajectory(self, coordinates, position=None, *, label, colors, dashed):
        trajectory = Trajectory(coordinates, position, label, colors, dashed)
        self._trajectories.append(trajectory)

        self._redraw_attractor()

        trace_coordinates, trace_position = self.__plot_coordinates_and_position(
            trajectory
        )

        return trace_coordinates, trace_position

    def _plot_trajectory(self, coordinates, *, label=None, color=None, trail=False):
        if self._attractor is None:
            raise ValueError(
                "An attractor must be set up first, please use "
                "set_attractor(Major_Body) or plot(orbit)"
            )

        colors = self._get_colors(color, trail)

        # Ensure that the coordinates are cartesian just in case,
        # to avoid weird errors later
        coordinates = coordinates.represent_as(CartesianRepresentation)

        return self.__add_trajectory(
            coordinates, None, label=str(label), colors=colors, dashed=False
        )

    def _plot(self, orbit, *, label=None, color=None, trail=False):
        colors = self._get_colors(color, trail)

        self.set_attractor(orbit.attractor)

        orbit = orbit.change_plane(self.plane)

        label = generate_label(orbit.epoch, label)
        coordinates = orbit.sample(self._num_points)

        return self.__add_trajectory(
            coordinates, orbit.r, label=label, colors=colors, dashed=True
        )

    def _plot_body_orbit(
        self,
        body,
        epoch,
        *,
        label=None,
        color=None,
        trail=False,
    ):
        if color is None:
            color = BODY_COLORS.get(body.name)

        self.set_attractor(body.parent)

        # Get approximate, mean value for the period
        period = get_mean_elements(body, epoch).period

        label = generate_label(epoch, label or str(body))
        epochs = time_range(
            epoch, periods=self._num_points, end=epoch + period, scale="tdb"
        )
        ephem = Ephem.from_body(body, epochs, attractor=body.parent, plane=self.plane)

        return self._plot_ephem(ephem, epoch, label=label, color=color, trail=trail)

    def _plot_ephem(self, ephem, epoch=None, *, label=None, color=None, trail=False):
        if self._attractor is None:
            raise ValueError(
                "An attractor must be set up first, please use "
                "set_attractor(Major_Body) or plot(orbit)"
            )

        if ephem.plane is not self.plane:
            raise ValueError(
                f"The ephemerides reference plane is {ephem.plane} "
                f"while the plotter is using {self.plane}, "
                "sample the ephemerides using a different plane "
                "or create a new plotter"
            )

        colors = self._get_colors(color, trail)

        coordinates = ephem.sample()
        if epoch is not None:
            r0 = ephem.rv(epoch)[0]
        else:
            r0 = None

        return self.__add_trajectory(
            coordinates, r0, label=str(label), colors=colors, dashed=False
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
        # Do not return the result of self._plot
        # This behavior might be overriden by subclasses
        self._plot_trajectory(coordinates, label=label, color=color, trail=trail)

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
        # Do not return the result of self._plot
        # This behavior might be overriden by subclasses
        self._plot(orbit, label=label, color=color, trail=trail)

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
        # Do not return the result of self._plot
        # This behavior might be overriden by subclasses
        self._plot_body_orbit(body, epoch, label=label, color=color, trail=trail)

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
        # Do not return the result of self._plot
        # This behavior might be overriden by subclasses
        self._plot_ephem(ephem, epoch, label=label, color=color, trail=trail)


class Mixin2D:
    _trajectories: List[Trajectory]

    def _redraw(self):
        raise NotImplementedError

    def _project(self, rr):
        rr_proj = rr - rr.dot(self._frame[2])[:, None] * self._frame[2]
        x = rr_proj.dot(self._frame[0])
        y = rr_proj.dot(self._frame[1])
        return x, y

    def _set_frame(self, p_vec, q_vec, w_vec):
        if not np.allclose([norm(v) for v in (p_vec, q_vec, w_vec)], 1):
            raise ValueError("Vectors must be unit.")
        elif not np.allclose([p_vec.dot(q_vec), q_vec.dot(w_vec), w_vec.dot(p_vec)], 0):
            raise ValueError("Vectors must be mutually orthogonal.")
        else:
            self._frame = p_vec, q_vec, w_vec

        if self._trajectories:
            self._redraw()

    def set_frame(self, p_vec, q_vec, w_vec):
        """Sets perifocal frame.

        Raises
        ------
        ValueError
            If the vectors are not a set of mutually orthogonal unit vectors.

        """
        warnings.warn(
            "Method set_frame is deprecated and will be removed in a future release, "
            "use `set_body_frame` or `set_orbit_frame` instead"
            "with your use case",
            DeprecationWarning,
            stacklevel=2,
        )
        self._set_frame(p_vec, q_vec, w_vec)

    def set_orbit_frame(self, orbit):
        """Sets perifocal frame based on an orbit.

        Parameters
        ----------
        orbit : ~poliastro.twobody.Orbit
            Orbit to use as frame.

        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            self._set_frame(*orbit.pqw())

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

        from poliastro.bodies import Sun
        from poliastro.twobody import Orbit

        from ..warnings import TimeScaleWarning

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
