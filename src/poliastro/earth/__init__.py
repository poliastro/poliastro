from typing import Dict

import numpy as np
from astropy import units as u

from poliastro.bodies import Earth
from poliastro.core.perturbations import J2_perturbation, atmospheric_drag_model
from poliastro.earth.enums import EarthGravity
from poliastro.twobody.orbit import Orbit
from poliastro.twobody.propagation import cowell


class EarthSatellite:

    """
    Position and velocity of a body with respect to Earth
    at a given time
    """

    def __init__(self, orbit):
        """Constructor.

        Parameters
        ----------
        orbit : Orbit
            Position and velocity of a body with respect to an attractor
            at a given time (epoch).

        Raises
        ------
        ValueError
            If the orbit's attractor is not Earth

        """

        if orbit.attractor is not Earth:
            raise ValueError("The attractor must be Earth")

        self._orbit = orbit  # type: Orbit

    @property
    def orbit(self):
        """Orbit of the EarthSatellite. """
        return self._orbit

    @classmethod
    @u.quantity_input(r=u.m, v=u.m / u.s)
    def from_vectors(cls, r, v):
        """Return EarthSatellite with the `Orbit` from position and velocity vectors.

        Parameters
        ----------
        r : ~astropy.units.Quantity
            Position vector wrt attractor center.
        v : ~astropy.units.Quantity
            Velocity vector.

        Returns
        -------
        EarthSatellite
            New EarthSatellite with the'Orbit' position and velocity vectors.


        """

        orbit = Orbit.from_vectors(Earth, r, v)
        return cls(orbit)

    @classmethod
    @u.quantity_input(a=u.m, ecc=u.one, inc=u.rad, raan=u.rad, argp=u.rad, nu=u.rad)
    def from_classical(
        cls, a, ecc, inc, raan, argp, nu,
    ):
        """Return EarthSatellite with the 'Orbit' from classical orbital elements.

        Parameters
        ----------
        a : ~astropy.units.Quantity
            Semi-major axis.
        ecc : ~astropy.units.Quantity
            Eccentricity.
        inc : ~astropy.units.Quantity
            Inclination
        raan : ~astropy.units.Quantity
            Right ascension of the ascending node.
        argp : ~astropy.units.Quantity
            Argument of the pericenter.
        nu : ~astropy.units.Quantity
            True anomaly.

        Returns
        -------
        EarthSatellite
            New EarthSatellite with the 'Orbit' from classical orbital elements.

        """

        orbit = Orbit.from_classical(Earth, a, ecc, inc, raan, argp, nu)
        return cls(orbit)

    def propagate(self, tof, atmosphere=None, gravity=None, A_over_m=None, *args):
        """Propagates an 'EarthSatellite Orbit' at a specified time.

        If value is true anomaly, propagate orbit to this anomaly and return the result.
        Otherwise, if time is provided, propagate this `EarthSatellite Orbit` some `time` and return the result.

        Parameters
        ----------

        tof : ~astropy.units.Quantity, ~astropy.time.Time, ~astropy.time.TimeDelta
            Scalar time to propagate.
        atmosphere:
            a callable model from poliastro.atmosphere
        gravity: EarthGravity
            There are two possible values, SPHERICAL and J2. Only J2 is implemented at the moment. Default value is None.
        A_over_m: float
            frontal area/mass of the spacecraft (km^2/kg)
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
        if atmosphere is not None and A_over_m is not None:
            perturbations[atmospheric_drag_model] = {
                "R": Earth.R.to(u.km).value,
                "C_D": 2.2,  # FIXME, add C_D as a parameter of the EarthSatellite object
                "A_over_m": A_over_m,
                "model": atmosphere,
            }
        ad_kwargs.update(perturbations=perturbations)
        new_orbit = self.orbit.propagate(value=tof, method=cowell, ad=ad, **ad_kwargs)
        return EarthSatellite(new_orbit)
