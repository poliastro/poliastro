from typing import Dict

import numpy as np
from astropy import units as u

from poliastro.bodies import Earth
from poliastro.constants import GM_earth
from poliastro.core.perturbations import J2_perturbation, atmospheric_drag_model
from poliastro.earth.enums import EarthGravity
from poliastro.spacecraft import Spacecraft
from poliastro.twobody.orbit import Orbit
from poliastro.twobody.propagation import cowell


class EarthSatellite:

    """
    Position and velocity of a body with respect to Earth
    at a given time
    """

    def __init__(self, orbit, spacecraft):
        """Constructor.

        Parameters
        ----------
        orbit : Orbit
            Position and velocity of a body with respect to an attractor
            at a given time (epoch).
        spacecraft: Spacecraft

        Raises
        ------
        ValueError
            If the orbit's attractor is not Earth

        """

        if orbit.attractor is not Earth:
            raise ValueError("The attractor must be Earth")

        self._orbit = orbit  # type: Orbit
        self._spacecraft = spacecraft  # type: Spacecraft

    @property
    def orbit(self):
        """Orbit of the EarthSatellite. """
        return self._orbit

    @property
    def spacecraft(self):
        """Spacecraft of the EarthSatellite. """
        return self._spacecraft

    @u.quantity_input(tof=u.min)
    def propagate(self, tof, atmosphere=None, gravity=None, *args):
        """Propagates an 'EarthSatellite Orbit' at a specified time.

        If value is true anomaly, propagate orbit to this anomaly and return the result.
        Otherwise, if time is provided, propagate this `EarthSatellite Orbit` some `time` and return the result.

        Parameters
        ----------

        tof : ~astropy.units.Quantity, ~astropy.time.Time, ~astropy.time.TimeDelta
            Scalar time to propagate.
        atmosphere:
            a callable model from poliastro.earth.atmosphere
        gravity: EarthGravity
            There are two possible values, SPHERICAL and J2. Only J2 is implemented at the moment. Default value is None.
        *args:
            parameters used in perturbation models.

        Returns
        -------
        EarthSatellite
            A new EarthSatellite with the propagated Orbit

        """

        ad_kwargs: Dict[object, dict] = {}
        perturbations: Dict[object, dict] = {}

        def ad(t0, state, k, perturbations):
            if perturbations:
                return np.sum(
                    [f(t0=t0, state=state, k=k, **p) for f, p in perturbations.items()],
                    axis=0,
                )
            else:
                return np.array([0, 0, 0])

        if gravity is EarthGravity.J2:
            perturbations[J2_perturbation] = {
                "J2": Earth.J2.value,
                "R": Earth.R.to(u.km).value,
            }
        if atmosphere is not None:
            perturbations[atmospheric_drag_model] = {
                "R": Earth.R.to(u.km).value,
                "C_D": self.spacecraft.C_D,
                "A_over_m": (self.spacecraft.A / self.spacecraft.m),
                "model": atmosphere,
            }
        ad_kwargs.update(perturbations=perturbations)
        new_orbit = self.orbit.propagate(value=tof, method=cowell, ad=ad, **ad_kwargs)
        return EarthSatellite(new_orbit, self.spacecraft)

    def rgt(self, ndays, norbits, atol=1e-10, iter=float("+inf")):
        """Calculates the orbit required for a repeating ground track orbit.

        Parameters
        ----------
        ndays: int
            Number of days.
        norbits: int
            Number of orbits in the repeat cycle.
        atol: float, optional
            Tolerance for the semi-major axis computation, default to 1e-8.
        iter: int
            Number of maximum desired iterations to obtain the most precise semi-major axis, default to infinite.

        Returns
        -------
        Orbit
            A new ground track orbit.

        Notes
        -----
        The algorithm was obtained from "Fundamentals of Astrodynamics and Applications", 4th ed (2013)" by David
        A. Vallado. For further information, please refer to pages from 875 to 879.
        This iterative numerical method is also described in “A Prograde Geosat Exact Repeat Mission?”,
        The Journal of the Astronautical Sciences, Vol. 39, No. 3, July-September 1991, pp. 313-326.

        """
        µ = GM_earth.to(u.km ** 3 / u.s ** 2)
        a = self.orbit.a
        e = self.orbit.ecc
        J2 = Earth.J2.value
        R = Earth.R.to(u.km)
        i = self.orbit.inc
        we = 7.292115855306589e-5 * (1 / u.s)

        ndays = int(ndays)
        norbits = int(norbits)
        k = norbits / ndays
        n = k * we

        atol = atol * u.km
        a_new = (µ / n ** 2) ** (1 / 3)
        e_new = e

        itt = 0
        while abs(a_new - a) >= atol and itt <= iter:
            a = a_new
            e = e_new
            p = a * (1 - e ** 2)
            factor = 1.5 * n * J2 * (R / p) ** 2
            dΩ = -factor * np.cos(i)
            dw = 0.5 * factor * (4 - 5 * np.sin(i) ** 2)
            dm = 0.5 * factor * (1 - e ** 2) ** 0.5 * (2 - 3 * np.sin(i) ** 2)
            n = k * (we - dΩ) - (dm + dw)
            a_new = (µ / n ** 2) ** (1 / 3)
            e_new = 1 - a / a_new * (1 - e)  # equivalent (new_a - rp)/new_a
            itt += 1

        new_orbit = Orbit.from_classical(
            Earth, a, e, i, self.orbit.raan, self.orbit.argp, self.orbit.nu
        )
        return new_orbit
