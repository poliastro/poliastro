from astropy import units as u
from matplotlib import pyplot as plt

from poliastro.plotting.util import generate_label


class GabbardPlotter:

    def __init__(self, ax, dark=False):

        self._ax = ax
        if not self._ax:
            if dark:
                with plt.style.context("dark_background"):
                    _, self._ax = plt.subplots(figsize=(6, 6))
            else:
                _, self._ax = plt.subplots(figsize=(6, 6))

        self._frame = None

    @staticmethod
    def _get_orbit_properties(orbit):
        return orbit.r_p, orbit.r_a, orbit.period

    def _get_orbit_property_list(self, orbits):
        apogees = []
        perigees = []
        periods = []

        for orbit in orbits:
            perigee, apogee, period = self._get_orbit_properties(orbit)
            apogees.append(apogee.value)
            perigees.append(perigee.value)
            periods.append(period.to(u.min).value)
        return apogees, perigees, periods

    def _static_gabbard_plot(self, orbits):
        """Plots a Static Gabbard Plot given a list of Orbits
        Parameters
        ----------
        orbits : ~poliastro.twobody.orbit.Orbit List
            The Orbits whose perigee and apogee will be plotted.
        """
        apogees, perigees, periods = self._get_orbit_property_list(orbits)

        apogee_paths = plt.scatter(periods, apogees, marker = "o", color = "blue", label = "Apogee")
        perigee_paths = plt.scatter(periods, perigees, marker = "o", color = "red", label = "Perigee")

        self._ax.set_xlabel("Period(min)")
        self._ax.set_ylabel("Altitude(km)")
        
        apogee_paths = plt.scatter(
            periods, apogees, marker="o", color="blue", label="Apogee"
        )
        perigee_paths = plt.scatter(
            periods, perigees, marker="o", color="red", label="Perigee"
        )

        self._ax.set_xlabel("Period(min)")
        self._ax.set_ylabel("Altitude(km)")

        return apogee_paths, perigee_paths

    def plot_orbits(self, orbits):
        apogee_paths, perigee_paths = self._static_gabbard_plot(orbits)
        self._set_legend(orbits[-1].epoch)
        return apogee_paths, perigee_paths

    def _set_legend(self, epoch):
        label = generate_label(epoch, "Epoch")
        print(label)
        if not self._ax.get_legend():
            size = self._ax.figure.get_size_inches() + [8, 0]
            self._ax.figure.set_size_inches(size)

        self._ax.legend(
            loc="upper left",
            bbox_to_anchor=(1.05, 1.015),
            title = label,
            numpoints=1,
        )
        
