""" Generates Tisserand plots """
from enum import Enum

import numpy as np
from astropy import units as u
from matplotlib import pyplot as plt

from poliastro.plotting._base import BODY_COLORS
from poliastro.twobody.mean_elements import get_mean_elements
from poliastro.util import norm


class TisserandKind(Enum):
    """All possible Tisserand kinds"""

    APSIS = "apsis"
    ENERGY = "energy"
    PERIOD = "period"


class TisserandPlotter:
    """Generates Tisserand figures"""

    def __init__(self, kind=TisserandKind.APSIS, axes=None):
        """Object initializer

        Parameters
        ----------
        kind : TisserandKind
            Nature for the Tisserand
        axes : ~matplotlib.pyplot.axes
            Axes for the figure

        """

        # Asign Tisserand kind
        self.kind = kind

        # Check if axis available
        if not axes:
            _, self.ax = plt.subplots(1, 1)
        else:
            self.ax = axes

        # Force axes scale regarding Tisserand kind
        self.ax.set_xscale("log")
        if self.kind == TisserandKind.APSIS:
            self.ax.set_yscale("log")

    def _solve_tisserand(
        self, body, vinf_span, num_contours, alpha_lim=(0, np.pi), N=100
    ):
        """Solves all possible Tisserand lines with a meshgrid workflow

        Parameters
        ----------
        body : ~poliastro.bodies.Body
            Body to be plotted Tisserand
        vinf_array : ~astropy.units.Quantity
            Desired Vinf for the flyby
        num_contours : int
            Number of contour lines for flyby speed
        alpha_lim : tuple
            Minimum and maximum flyby angles.
        N : int
            Number of points for flyby angle.

        Notes
        -----
        The algorithm for generating Tisserand plots is the one depicted in
        "Preliminary Trajectory Design of a Mission to Enceladus" by David
        Falcato Fialho Palma, section 3.6

        """

        # Generate mean orbital elements Earth
        body_rv = get_mean_elements(body).to_vectors()
        R_body, V_body = norm(body_rv.r), norm(body_rv.v)

        # Generate non-dimensional velocity and alpha span
        vinf_array = np.linspace(vinf_span[0], vinf_span[-1], num_contours)
        alpha_array = np.linspace(alpha_lim[0], alpha_lim[-1], N)
        vinf_array /= V_body

        # Construct the mesh for any configuration
        V_INF, ALPHA = np.meshgrid(vinf_array, alpha_array)

        # Solving for non-dimensional a_sc and ecc_sc
        A_SC = 1 / np.abs(1 - V_INF**2 - 2 * V_INF * np.cos(ALPHA))
        ECC_SC = np.sqrt(
            1 - 1 / A_SC * ((3 - 1 / A_SC - V_INF**2) / (2)) ** 2
        )

        # Compute main Tisserand variables
        RR_P = A_SC * R_body * (1 - ECC_SC)
        RR_A = A_SC * R_body * (1 + ECC_SC)
        TT = 2 * np.pi * np.sqrt((A_SC * R_body) ** 3 / body.parent.k)
        EE = -body.parent.k / (2 * A_SC * R_body)

        # Build color lines to internal canvas
        return RR_P, RR_A, EE, TT

    def _build_lines(self, RR_P, RR_A, EE, TT, color):
        """Collect lines and append them to internal data

        Parameters
        ----------
        data : list
            Array containing [RR_P, RR_A, EE, TT, color]

        Returns
        -------
        lines: list
            Plotting lines for the Tisserand
        """

        # Plot desired kind lines
        if self.kind == TisserandKind.APSIS:
            # Generate apsis lines
            lines = self.ax.plot(RR_A.to(u.AU), RR_P.to(u.AU), color=color)
        elif self.kind == TisserandKind.ENERGY:
            # Generate energy lines
            lines = self.ax.plot(
                RR_P.to(u.AU), EE.to(u.km**2 / u.s**2), color=color
            )
        elif self.kind == TisserandKind.PERIOD:
            # Generate period lines
            lines = self.ax.plot(RR_P.to(u.AU), TT.to(u.year), color=color)

        return lines

    def plot_line(self, body, vinf, alpha_lim=(0, np.pi), color=None):
        """Plots body Tisserand line within flyby angle

        Parameters
        ----------
        body : ~poliastro.bodies.Body
            Body to be plotted Tisserand
        vinf : ~astropy.units.Quantity
            Vinf velocity line
        alpha_lim : tuple
            Minimum and maximum flyby angles
        color : str
            String representing for the color lines

        Returns
        -------
        self.ax: ~matplotlib.axes.Axes
            Apsis tisserand is the default plotting option

        """

        # HACK: to reuse Tisserand solver, we transform input Vinf into a tuple
        vinf_span = (vinf, vinf)

        # Solve Tisserand parameters
        RR_P, RR_A, EE, TT = self._solve_tisserand(
            body, vinf_span, num_contours=2, alpha_lim=alpha_lim
        )

        # Check if color defined
        if not color:
            color = BODY_COLORS[body.name]

        # Build canvas lines from Tisserand parameters
        self._build_lines(RR_P, RR_A, EE, TT, color)

        return self.ax

    def plot(self, body, vinf_span, num_contours=10, color=None):
        """Plots body Tisserand for given amount of solutions within Vinf span

        Parameters
        ----------
        body : ~poliastro.bodies.Body
            Body to be plotted Tisserand
        vinf_span : tuple
            Minimum and maximum Vinf velocities
        num_contours : int
            Number of points to iterate over previously defined velocities
        color : str
            String representing for the color lines

        Returns
        -------
        self.ax: ~matplotlib.axes.Axes
            Apsis tisserand is the default plotting option

        """

        # Solve Tisserand parameters
        RR_P, RR_A, EE, TT = self._solve_tisserand(
            body, vinf_span, num_contours
        )

        # Check if color defined
        if not color:
            color = BODY_COLORS[body.name]

        # Build canvas lines from Tisserand parameters
        self._build_lines(RR_P, RR_A, EE, TT, color)

        return self.ax
