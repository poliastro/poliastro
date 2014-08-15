# coding: utf-8
""" Plotting utilities.

Notes
-----
TODO: I still miss a way to plot several orbits in one plot.

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
    def __init__(self, num_points=100):
        """Constructor.

        """
        self.num_points = num_points

    def plot(self, state, ax, osculating=True):
        """Plots state and osculating orbit in their plane.

        """
        if not ax:
            _, ax = plt.subplots(figsize=(6, 6))

        lines = []
        num = self.num_points
        # FIXME Some faulty logic here
        nu_vals = np.linspace(0, 2 * np.pi, num) + state.nu.to(u.rad).value
        r_pqw, _ = rv_pqw(state.attractor.k.to(u.km ** 3 / u.s ** 2).value,
                          state.p.to(u.km).value, state.ecc.value, nu_vals)

        # Current position
        l, = ax.plot(r_pqw[0, 0], r_pqw[0, 1], 'o')
        lines.append(l)

        # Attractor
        attractor = mpl.patches.Circle((0, 0),
                                       state.attractor.R.to(u.km).value,
                                       lw=0, color='#204a87')  # Earth
        ax.add_patch(attractor)

        if osculating:
            l, = ax.plot(r_pqw[:, 0], r_pqw[:, 1], '--', color=l.get_color())
            lines.append(l)

        # Appearance
        ax.set_title(state.epoch.iso)
        ax.set_xlabel("$x$ (km)")
        ax.set_ylabel("$y$ (km)")
        ax.set_aspect(1)

        return lines
