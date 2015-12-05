# coding: utf-8
""" Plotting utilities.

"""

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

from astropy import units as u

from poliastro.twobody.classical import rv_pqw
from poliastro.util import norm


def plot(state):
    """Plots a ``State``.

    For more advanced tuning, use the :py:class:`OrbitPlotter` class.

    """
    op = OrbitPlotter()
    return op.plot(state)


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
        ax : Axes
            Axes in which to plot. If not given, new ones will be created.
        num_points : int, optional
            Number of points to use in plots, default to 100.

        """
        self.ax = ax
        if not self.ax:
            _, self.ax = plt.subplots(figsize=(6, 6))
        self.num_points = num_points
        self._frame = None
        self._states = []

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

    def plot(self, state, osculating=True):
        """Plots state and osculating orbit in their plane.

        """
        if not self._frame:
            self.set_frame(*state.pqw())

        self._states.append(state)

        lines = []

        # Generate points of the osculating orbit
        if state.ecc >= 1.0:
            # Select a sensible limiting value for non-closed orbits.
            # This corresponds to r = 3p.
            max_nu = np.arccos(-(1 - 1 / 3.) / state.ecc).to(u.rad).value
        else:
            max_nu = np.pi  # rad
        nu_vals = np.linspace(-max_nu, max_nu, self.num_points)  # rad

        # Insert state true anomaly into array
        idx = np.searchsorted(nu_vals, state.nu.to(u.rad))
        nu_vals = np.insert(nu_vals, idx,
                            state.nu.to(u.rad).value) * u.rad

        # Compute PQW coordinates
        r_pqw, _ = rv_pqw(state.attractor.k.to(u.km ** 3 / u.s ** 2).value,
                          state.p.to(u.km).value,
                          state.ecc.value,
                          nu_vals.value)
        r_pqw = r_pqw * u.km

        # Express on inertial frame
        e_vec, p_vec, h_vec = state.pqw()
        p_vec = np.cross(h_vec, e_vec) * u.one
        rr = (r_pqw[:, 0, None].dot(e_vec[None, :]) +
              r_pqw[:, 1, None].dot(p_vec[None, :]))

        # Project on OrbitPlotter frame
        # x_vec, y_vec, z_vec = self._frame
        rr_proj = rr - rr.dot(self._frame[2])[:, None] * self._frame[2]
        x = rr_proj.dot(self._frame[0])
        y = rr_proj.dot(self._frame[1])

        # Plot current position
        l, = self.ax.plot(x[idx].to(u.km).value, y[idx].to(u.km).value,
                          'o', mew=0)
        lines.append(l)

        # Attractor
        self.ax.add_patch(mpl.patches.Circle((0, 0),
                                             state.attractor.R.to(u.km).value,
                                             lw=0, color='#204a87'))  # Earth

        if osculating:
            l, = self.ax.plot(x.to(u.km).value, y.to(u.km).value,
                              '--', color=l.get_color())
            lines.append(l)

        self.ax.set_title(state.epoch.iso)
        self.ax.set_xlabel("$x$ (km)")
        self.ax.set_ylabel("$y$ (km)")
        self.ax.set_aspect(1)

        return lines
