""" Plotting utilities.

"""
from typing import List

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

from plotly.offline import iplot
from plotly.graph_objs import Scatter3d, Surface, Layout

from astropy import units as u
from astropy.coordinates.angles import Angle

from poliastro.twobody.classical import rv_pqw
from poliastro.util import norm


BODY_COLORS = {
    "Sun": "#ffcc00",
    "Earth": "#204a87",
    "Jupiter": "#ba3821",
}


def plot(state, label=None, color=None):
    """Plots a ``State``.

    For more advanced tuning, use the :py:class:`OrbitPlotter` class.

    """
    op = OrbitPlotter()
    return op.plot(state, label=label, color=color)


def plot3d(orbit, sampling=None, *, label=None, color=None):
    frame = OrbitPlotter3D()
    frame.plot(orbit, sampling, label=label, color=color)
    frame.show()

    return frame


class OrbitPlotter(object):
    """OrbitPlotter class.

    This class holds the perifocal plane of the first
    :py:class:`~poliastro.twobody.State` plotted in it using
    :py:meth:`plot`, so all following
    plots will be projected on that plane. Alternatively, you can call
    :py:meth:`set_frame` to set the frame before plotting.

    """
    def __init__(self, ax=None, num_points=100):
        """Constructor.

        Parameters
        ----------
        ax : ~matplotlib.axes.Axes
            Axes in which to plot. If not given, new ones will be created.
        num_points : int, optional
            Number of points to use in plots, default to 100.

        """
        self.ax = ax
        if not self.ax:
            _, self.ax = plt.subplots(figsize=(6, 6))
        self.num_points = num_points
        self._frame = None
        self._attractor_radius = None

    def set_frame(self, p_vec, q_vec, w_vec):
        """Sets perifocal frame if not existing.

        Raises
        ------
        ValueError
            If the vectors are not a set of mutually orthogonal unit vectors.
        NotImplementedError
            If the frame was already set up.

        """
        if not self._frame:
            if not np.allclose([norm(v) for v in (p_vec, q_vec, w_vec)], 1):
                raise ValueError("Vectors must be unit.")
            if not np.allclose([p_vec.dot(q_vec),
                                q_vec.dot(w_vec),
                                w_vec.dot(p_vec)], 0):
                raise ValueError("Vectors must be mutually orthogonal.")
            else:
                self._frame = p_vec, q_vec, w_vec
        else:
            raise NotImplementedError

    def set_attractor(self, orbit):
        """Sets plotting attractor.

        Parameters
        ----------
        orbit : ~poliastro.twobody.orbit.Orbit
            orbit with attractor to plot.

        """
        radius = max(orbit.attractor.R.to(u.km).value,
                     orbit.r_p.to(u.km).value / 6)
        color = BODY_COLORS.get(orbit.attractor.name, "#999999")
        self._attractor_radius = radius

        for attractor in self.ax.findobj(match=mpl.patches.Circle):
            attractor.remove()

        self.ax.add_patch(mpl.patches.Circle((0, 0), radius, lw=0, color=color))

    def plot(self, orbit, label=None, color=None):
        """Plots state and osculating orbit in their plane.

        """
        # TODO: This function needs a refactoring
        if not self._frame:
            self.set_frame(*orbit.pqw())

        # if new attractor radius is smaller, plot it
        new_radius = max(orbit.attractor.R.to(u.km).value,
                         orbit.r_p.to(u.km).value / 6)
        if not self._attractor_radius:
            self.set_attractor(orbit)
        elif new_radius < self._attractor_radius:
            self.set_attractor(orbit)

        lines = []

        nu_vals = self._generate_vals(orbit.state)

        # Compute PQW coordinates
        r_pqw, _ = rv_pqw(orbit.attractor.k.to(u.km ** 3 / u.s ** 2).value,
                          orbit.p.to(u.km).value,
                          orbit.ecc.value,
                          nu_vals.value)
        r_pqw = r_pqw * u.km

        # Express on inertial frame
        e_vec, p_vec, h_vec = orbit.pqw()
        p_vec = np.cross(h_vec, e_vec) * u.one
        rr = (r_pqw[:, 0, None].dot(e_vec[None, :]) +
              r_pqw[:, 1, None].dot(p_vec[None, :]))

        # Project on OrbitPlotter frame
        # x_vec, y_vec, z_vec = self._frame
        rr_proj = rr - rr.dot(self._frame[2])[:, None] * self._frame[2]
        x = rr_proj.dot(self._frame[0])
        y = rr_proj.dot(self._frame[1])

        # Plot current position
        l, = self.ax.plot(x[0].to(u.km).value, y[0].to(u.km).value,
                          'o', mew=0, color=color)
        lines.append(l)

        l, = self.ax.plot(x.to(u.km).value, y.to(u.km).value,
                          '--', color=l.get_color())
        lines.append(l)

        if label:
            # This will apply the label to either the point or the osculating
            # orbit depending on the last plotted line, as they share variable
            if not self.ax.get_legend():
                size = self.ax.figure.get_size_inches() + [8, 0]
                self.ax.figure.set_size_inches(size)
            orbit.epoch.out_subfmt = 'date_hm'
            label = '{} ({})'.format(orbit.epoch.iso, label)

            l.set_label(label)
            self.ax.legend(bbox_to_anchor=(1.05, 1), title="Names and epochs")

        self.ax.set_xlabel("$x$ (km)")
        self.ax.set_ylabel("$y$ (km)")
        self.ax.set_aspect(1)

        return lines

    def _generate_vals(self, state):
        """Generate points of the osculating orbit.

        """
        nu_vals = Angle((np.linspace(0, 2 * np.pi, self.num_points) +
                         state.nu.to(u.rad).value) * u.rad).wrap_at('360d')

        if state.ecc >= 1:
            # Select a sensible limiting value for non-closed orbits
            # This corresponds to r = 3p
            nu_limit = Angle(np.arccos(-(1 - 1 / 3.) / state.ecc))
            nu_invalid = ((nu_vals > nu_limit) &
                          (nu_vals < (-nu_limit).wrap_at('360d')))
            nu_vals[nu_invalid] = np.nan

        return nu_vals


def _generate_label(orbit, label):
    orbit.epoch.out_subfmt = 'date_hm'
    label_ = '{}'.format(orbit.epoch.iso)
    if label:
        orbit.epoch.out_subfmt = 'date_hm'
        label_ += ' ({})'.format(label)

    return label_


def _generate_sphere(radius, num=20):
    u1 = np.linspace(0, 2 * np.pi, num)
    v1 = u1.copy()
    uu, vv = np.meshgrid(u1, v1)

    xx = radius * np.cos(uu) * np.sin(vv)
    yy = radius * np.sin(uu) * np.sin(vv)
    zz = radius * np.cos(vv)

    return xx, yy, zz


class OrbitPlotter3D:
    def __init__(self):
        self._orbits = []  # type: List[tuple]
        self._attractor_radius = np.inf

        self._layout = Layout(
            autosize=False,
            scene=dict(
                xaxis=dict(
                    title="x (km)",
                ),
                yaxis=dict(
                    title="y (km)",
                ),
                zaxis=dict(
                    title="z (km)",
                ),
                aspectmode="data",  # Important!
            ),
        )
        self._data = []  # type: List[dict]

    def _draw_attractor(self, orbit, threshold=0.15):
        # Select a sensible value for the radius: realistic for low orbits,
        # visible for high and very high orbits
        radius = max(orbit.attractor.R.to(u.km), orbit.r_p.to(u.km) * threshold).value

        # If the resulting radius is smaller than the current one, redraw it
        if radius < self._attractor_radius:
            name = orbit.attractor.name
            color = BODY_COLORS.get(name, "#999999")

            xx, yy, zz = _generate_sphere(radius)

            sphere = Surface(
                x=xx, y=yy, z=zz,
                name=name,
                colorscale=[[0, color], [1, color]],
                cauto=False, cmin=1, cmax=1, showscale=False,  # Boilerplate
            )

            self._attractor_radius = radius
            self._data.append(sphere)

    def _plot_orbit(self, orbit, sampling, label, color):
        rr = orbit.sample(sampling)
        trace = Scatter3d(
            x=rr.x.value, y=rr.y.value, z=rr.z.value,
            name=label,
            line=dict(
                color=color,
                width=5,
            ),
            mode="lines",  # Boilerplate
        )
        self._data.append(trace)

    def plot(self, orbit, sampling=None, *, label=None, color=None):
        label = _generate_label(orbit, label)

        self._orbits.append(
            (orbit, sampling, label, color)
        )

        self._draw_attractor(orbit)

        self._plot_orbit(orbit, sampling, label, color)

    def show(self):
        iplot(dict(
            data=self._data,
            layout=self._layout,
        ))
