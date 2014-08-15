# coding: utf-8
""" Plotting utilities.

"""

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

from astropy import units as u
u.one = u.dimensionless_unscaled  # astropy #1980

from poliastro.twobody.conversion import rv_pqw


class OrbitPlotter(object):
    """OrbitPlotter class.

    """
    def __init__(self, ax, num_points=100):
        """Constructor.

        """
        self.ax = ax
        if not self.ax:
            _, self.ax = plt.subplots(figsize=(6, 6))
        self.num_points = num_points
        self._frame = None
        self._states = []

    def plot(self, state, osculating=True):
        """Plots state and osculating orbit in their plane.

        """
        if not self._frame:
            self._frame = state.pqw()
        self._states.append(state)

        lines = []

        # Generate points of the osculating orbit
        nu_vals = (np.linspace(0, 2 * np.pi, self.num_points) +
                   state.nu.to(u.rad).value)
        r_pqw, _ = rv_pqw(state.attractor.k.to(u.km ** 3 / u.s ** 2).value,
                          state.p.to(u.km).value, state.ecc.value, nu_vals)

        # Express on inertial frame
        e_vec, p_vec, h_vec = state.pqw()
        p_vec = np.cross(h_vec, e_vec) * u.one
        rr = (r_pqw[:, 0, None].dot(e_vec[None, :]) +
              r_pqw[:, 1, None].dot(p_vec[None, :])) * u.km

        # Project on OrbitPlotter frame
        x_vec, y_vec, z_vec = self._frame
        rr_proj = rr - rr.dot(z_vec)[:, None] * z_vec
        x = rr_proj.dot(x_vec)
        y = rr_proj.dot(y_vec)

        # Plot current position
        l, = self.ax.plot(x[0].to(u.km).value, y[0].to(u.km).value,
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
