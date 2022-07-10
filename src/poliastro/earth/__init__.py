""" Earth focused orbital mechanics routines """

from typing import Dict

import numpy as np
from astropy import units as u

from poliastro.bodies import Earth
from poliastro.core.perturbations import J2_perturbation
from poliastro.core.propagation import func_twobody
from poliastro.earth.enums import EarthGravity
from poliastro.spacecraft import Spacecraft
from poliastro.twobody.orbit import Orbit
from poliastro.twobody.propagation import CowellPropagator


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
        spacecraft : Spacecraft

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
        """Orbit of the EarthSatellite."""
        return self._orbit

    @property
    def spacecraft(self):
        """Spacecraft of the EarthSatellite."""
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
        gravity : EarthGravity
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
                    [
                        f(t0=t0, state=state, k=k, **p)
                        for f, p in perturbations.items()
                    ],
                    axis=0,
                )
            else:
                return np.array([0, 0, 0])

        if gravity is EarthGravity.J2:
            perturbations[J2_perturbation] = {
                "J2": Earth.J2.value,
                "R": Earth.R.to_value(u.km),
            }
        if atmosphere is not None:
            # Cannot compute density without knowing the state,
            # the perturbations parameters are not always fixed
            # TODO: This whole function probably needs a refactoring
            raise NotImplementedError

        def f(t0, state, k):
            du_kep = func_twobody(t0, state, k)
            ax, ay, az = ad(t0, state, k, perturbations)
            du_ad = np.array([0, 0, 0, ax, ay, az])

            return du_kep + du_ad

        ad_kwargs.update(perturbations=perturbations)
        new_orbit = self.orbit.propagate(tof, method=CowellPropagator(f=f))
        return EarthSatellite(new_orbit, self.spacecraft)
