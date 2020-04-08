import warnings
from collections import namedtuple
from typing import List

import numpy as np
from astropy import units as u
from astropy.coordinates import CartesianRepresentation

from poliastro.plotting.util import BODY_COLORS, generate_label
from poliastro.util import norm

from ..frames import Planes


class Trajectory(namedtuple("Trajectory", ["positions", "state", "label", "color"])):
    pass


class BaseOrbitPlotter:
    """
    Base class for all the OrbitPlotter classes.
    """

    def __init__(self, num_points=150):
        self._num_points = num_points

        self._trajectories = []  # type: List[Trajectory]

        self._attractor = None
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

    def _set_plane(self, plane, fail_if_set=True):
        if self._plane is None:
            self._plane = plane
        elif plane is not self._plane and fail_if_set:
            raise NotImplementedError(f"Plane has already been set to {self._plane}")

    def set_attractor(self, attractor):
        """Sets plotting attractor.

        Parameters
        ----------
        attractor : ~poliastro.bodies.Body
            Central body.

        """
        self._set_attractor(attractor)

    def set_plane(self, plane):
        """Sets reference plane.

        Parameters
        ----------
        plane : ~poliastro.frames.enums.Planes
            Reference plane.

        """
        self._set_plane(plane)

    def _clear_attractor(self):
        raise NotImplementedError

    def _redraw_attractor(self):
        # Select a sensible value for the radius: realistic for low orbits,
        # visible for high and very high orbits
        min_radius = min(
            [
                positions.represent_as(CartesianRepresentation).norm().min() * 0.15
                for positions, _, _, _ in self._trajectories
            ]
            or [0 * u.m]
        )
        radius = max(self._attractor.R.to(u.km), min_radius.to(u.km))

        color = BODY_COLORS.get(self._attractor.name, "#999999")

        self._clear_attractor()

        if radius < self._attractor_radius:
            self._attractor_radius = radius

        self._draw_sphere(
            self._attractor_radius, color, self._attractor.name,
        )

    def _get_colors(self, color, trail):
        raise NotImplementedError

    def _draw_point(self, radius, color, name, center=None):
        raise NotImplementedError

    def _draw_sphere(self, radius, color, name, center=None):
        raise NotImplementedError

    def _plot_trajectory(self, positions, label, colors, dashed):
        raise NotImplementedError

    def _plot_r(self, state, label, colors):
        radius = min(
            self._attractor_radius * 0.5, (norm(state) - self._attractor.R) * 0.5
        )  # Arbitrary thresholds
        self._draw_point(radius, colors[0], label, center=state)

    def _plot(self, positions, state, label, colors):
        trace_trajectory = self._plot_trajectory(positions, label, colors, True)

        # Redraw the attractor now to compute the attractor radius
        # with the trajectory we just added
        self._redraw_attractor()

        if state is not None:
            # Plot required 2D/3D shape in the position of the body
            trace_r = self._plot_r(state, label, colors)
        else:
            trace_r = None

        return trace_trajectory, trace_r

    def plot_trajectory(self, positions, *, label=None, color=None, trail=False):
        """Plots a precomputed trajectory.

        An attractor must be set first.

        Parameters
        ----------
        positions : ~astropy.coordinates.CartesianRepresentation
            Trajectory to plot.
        label : string, optional
            Label of the trajectory.
        color : string, optional
            Color of the trajectory.
        trail : bool, optional
            Fade the orbit trail, default to False.

        """
        if self._attractor is None:
            raise ValueError(
                "An attractor must be set up first, please use "
                "set_attractor(Major_Body) or plot(orbit)"
            )

        colors = self._get_colors(color, trail)

        self._plot(positions, None, str(label), colors)

        self._trajectories.append(Trajectory(positions, None, str(label), colors[0]))

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
        colors = self._get_colors(color, trail)

        self.set_attractor(orbit.attractor)
        # If plane is already set, we will use the current one to reproject
        self._set_plane(orbit.plane, fail_if_set=False)

        label = generate_label(orbit.epoch, label)
        positions = orbit.change_plane(self._plane).sample(self._num_points)

        self._trajectories.append(Trajectory(positions, orbit.r, label, colors[0]))

        self._plot(positions, orbit.r, label, colors)

    def plot_body_orbit(
        self,
        body,
        epoch=None,
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
        epoch : astropy.time.Time, optional
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
        from poliastro.twobody import Orbit

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            orbit = Orbit.from_body_ephem(body, epoch)

        self.plot(orbit, label=label or str(body), color=color, trail=trail)


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
        self._set_frame(*orbit.pqw())

    def set_body_frame(self, body, epoch=None):
        """Sets perifocal frame based on the orbit of a body at a particular epoch if given.

        Parameters
        ----------
        body : poliastro.bodies.SolarSystemBody
            Body.
        epoch : astropy.time.Time, optional
            Epoch of current position.

        """
        from poliastro.twobody import Orbit

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            orbit = Orbit.from_body_ephem(body, epoch)

        self.set_orbit_frame(orbit)
