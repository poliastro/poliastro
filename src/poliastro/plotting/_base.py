from collections import namedtuple
from typing import List

import numpy as np
from astropy import units as u
from astropy.coordinates import CartesianRepresentation

from poliastro.plotting.util import BODY_COLORS, generate_label
from poliastro.util import norm


class Trajectory(namedtuple("Trajectory", ["positions", "state", "label", "color"])):
    pass


class BaseOrbitPlotter:
    """
    Parent Class for the 2D and 3D OrbitPlotter Classes based on Plotly.
    """

    def __init__(self, num_points=150):
        self._num_points = num_points

        self._trajectories = []  # type: List[Trajectory]

        self._attractor = None

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

        trace = self._plot_trajectory(positions, str(label), colors, False)

        self._trajectories.append(Trajectory(positions, None, label, colors[0]))

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

        self._set_attractor(orbit.attractor)

        label = generate_label(orbit, label)
        positions = orbit.sample(self._num_points)

        trace = self._plot_trajectory(positions, label, colors, True)

        self._trajectories.append(Trajectory(positions, orbit.r, label, colors[0]))

        # Redraw the attractor now to compute the attractor radius
        self._redraw_attractor()

        # Plot required 2D/3D shape in the position of the body
        radius = min(
            self._attractor_radius * 0.5, (norm(orbit.r) - orbit.attractor.R) * 0.5
        )  # Arbitrary thresholds
        self._draw_point(radius, colors[0], label, center=orbit.r)


class Mixin2D:
    def _project(self, rr):
        rr_proj = rr - rr.dot(self._frame[2])[:, None] * self._frame[2]
        x = rr_proj.dot(self._frame[0])
        y = rr_proj.dot(self._frame[1])
        return x, y

    def set_frame(self, p_vec, q_vec, w_vec):
        """Sets perifocal frame.

        Raises
        ------
        ValueError
            If the vectors are not a set of mutually orthogonal unit vectors.

        """
        if not np.allclose([norm(v) for v in (p_vec, q_vec, w_vec)], 1):
            raise ValueError("Vectors must be unit.")
        elif not np.allclose([p_vec.dot(q_vec), q_vec.dot(w_vec), w_vec.dot(p_vec)], 0):
            raise ValueError("Vectors must be mutually orthogonal.")
        else:
            self._frame = p_vec, q_vec, w_vec

        if self._trajectories:
            self._redraw()
