import warnings
from collections import namedtuple
from typing import List

import numpy as np
from astropy import units as u
from astropy.coordinates import CartesianRepresentation

from poliastro.plotting.util import BODY_COLORS, generate_label
from poliastro.util import norm

from ..frames import Planes


class Trajectory(
    namedtuple("Trajectory", ["coordinates", "position", "label", "color"])
):
    pass


class BaseOrbitPlotter:
    """
    Base class for all the OrbitPlotter classes.
    """

    def __init__(self, num_points=150):
        self._num_points = num_points

        self._trajectories = []  # type: List[Trajectory]

        self._attractor = None

        # This plane is used as a reference
        # to conceal orbits in different planes,
        # but it is not exposed as public API
        # because it can be confusing with the 2D frames
        self._plane = None

        self._attractor_radius = np.inf * u.km

    @property
    def trajectories(self):
        return self._trajectories

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
            [coordinates.norm().min() for coordinates, _, _, _ in self._trajectories]
            or [0 * u.m]
        )
        self._attractor_radius = max(
            self._attractor.R.to(u.km), min_distance.to(u.km) * 0.15
        )

        color = BODY_COLORS.get(self._attractor.name, "#999999")

        self._clear_attractor()

        self._draw_sphere(
            self._attractor_radius, color, self._attractor.name,
        )

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

        self._trajectories.append(Trajectory(coordinates, None, str(label), colors))

        self._redraw_attractor()

        trace_coordinates = self._plot_coordinates(coordinates, label, colors, False)

        return trace_coordinates

    def _plot(self, orbit, *, label=None, color=None, trail=False):
        colors = self._get_colors(color, trail)

        self.set_attractor(orbit.attractor)

        if self._plane is None:
            self._plane = orbit.plane
        elif orbit.plane is not self._plane:
            orbit = orbit.change_plane(self._plane)

        label = generate_label(orbit.epoch, label)
        coordinates = orbit.sample(self._num_points)

        self._trajectories.append(Trajectory(coordinates, orbit.r, str(label), colors))

        self._redraw_attractor()

        trace_coordinates = self._plot_coordinates(coordinates, label, colors, True)
        trace_position = self._plot_position(orbit.r, label, colors)

        return trace_coordinates, trace_position

    def _plot_body_orbit(
        self,
        body,
        epoch,
        plane=Planes.EARTH_ECLIPTIC,
        *,
        label=None,
        color=None,
        trail=False,
    ):
        if self._plane is None:
            self._plane = plane

        if color is None:
            color = BODY_COLORS.get(body.name)

        from poliastro.twobody import Orbit

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            orbit = Orbit.from_body_ephem(body, epoch)

        return self._plot(orbit, label=label or str(body), color=color, trail=trail)

    def plot_trajectory(self, coordinates, *, label=None, color=None, trail=False):
        """Plots a precomputed trajectory.

        An attractor must be set first.

        Parameters
        ----------
        coordinates : ~astropy.coordinates.CartesianRepresentation
            Trajectory to plot.
        label : string, optional
            Label of the trajectory.
        color : string, optional
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
        label : string, optional
            Label of the orbit.
        color : string, optional
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
        plane=Planes.EARTH_ECLIPTIC,
        *,
        label=None,
        color=None,
        trail=False,
    ):
        """Plots complete revolution of body and current position.

        Parameters
        ----------
        body : poliastro.bodies.SolarSystemBody
            Body.
        epoch : astropy.time.Time
            Epoch of current position.
        plane : ~poliastro.frames.enums.Planes
            Reference plane.
        label : str, optional
            Label of the orbit, default to the name of the body.
        color : string, optional
            Color of the line and the position.
        trail : bool, optional
            Fade the orbit trail, default to False.

        """
        # Do not return the result of self._plot
        # This behavior might be overriden by subclasses
        self._plot_body_orbit(body, epoch, plane, label=label, color=color, trail=trail)


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

    def set_body_frame(self, body, epoch=None, plane=Planes.EARTH_ECLIPTIC):
        """Sets perifocal frame based on the orbit of a body at a particular epoch if given.

        Parameters
        ----------
        body : poliastro.bodies.SolarSystemBody
            Body.
        epoch : astropy.time.Time, optional
            Epoch of current position.
        plane : ~poliastro.frames.enums.Planes
            Reference plane.

        """
        from poliastro.twobody import Orbit

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            orbit = Orbit.from_body_ephem(body, epoch).change_plane(plane)

        self.set_orbit_frame(orbit)
