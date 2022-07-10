from astropy import units as u
from matplotlib import pyplot as plt

from poliastro.plotting.util import generate_label


class GabbardPlotter:
    """GabbardPlotter class."""

    def __init__(
        self, ax=None, dark=False, altitude_unit=u.km, period_unit=u.min
    ):

        self._ax = ax
        if not self._ax:
            if dark:
                with plt.style.context("dark_background"):
                    _, self._ax = plt.subplots(figsize=(6, 6))
            else:
                _, self._ax = plt.subplots(figsize=(6, 6))

        self._frame = None
        self._altitude_unit = altitude_unit
        self._period_unit = period_unit

    def _get_orbit_property_list(self, orbits):
        apogees = []
        perigees = []
        periods = []

        for orbit in orbits:
            perigee, apogee, period = orbit.r_p, orbit.r_a, orbit.period
            apogees.append(apogee.to_value(self._altitude_unit))
            perigees.append(perigee.to_value(self._altitude_unit))
            periods.append(period.to_value(self._period_unit))
        return apogees, perigees, periods

    def _static_gabbard_plot(self, orbits):
        """Plots a Static Gabbard Plot given a list of Orbits

        Parameters
        ----------
        orbits : ~poliastro.twobody.orbit.Orbit List
            The Orbits whose perigee and apogee will be plotted.

        """
        apogees, perigees, periods = self._get_orbit_property_list(orbits)

        apogee_paths = plt.scatter(
            periods, apogees, marker="o", color="blue", label="Apogee"
        )
        perigee_paths = plt.scatter(
            periods, perigees, marker="o", color="red", label="Perigee"
        )

        self._ax.set_xlabel(f"Period ({self._period_unit:s})")
        self._ax.set_ylabel(f"Altitude ({self._altitude_unit:s})")

        return apogee_paths, perigee_paths

    def plot_orbits(self, orbits, label=""):
        apogee_paths, perigee_paths = self._static_gabbard_plot(orbits)
        self._set_legend(orbits[-1].epoch, label)
        return apogee_paths, perigee_paths

    def _set_legend(self, epoch, label):
        label = generate_label(epoch, label)
        if not self._ax.get_legend():
            size = self._ax.figure.get_size_inches() + [8, 0]
            self._ax.figure.set_size_inches(size)

        self._ax.legend(
            loc="upper left",
            bbox_to_anchor=(1.05, 1.015),
            title=label,
            numpoints=1,
        )
