from warnings import warn

import numpy as np

from astropy import units as u
from astropy import time
from astropy.coordinates import CartesianRepresentation, get_body_barycentric_posvel

from poliastro.constants import J2000
from poliastro.twobody.angles import nu_to_M
from poliastro.twobody.propagation import propagate, cowell

from poliastro.twobody import rv
from poliastro.twobody import classical
from poliastro.twobody import equinoctial

from ._base import BaseState


ORBIT_FORMAT = "{r_p:.0f} x {r_a:.0f} x {inc:.1f} orbit around {body}"


class TimeScaleWarning(UserWarning):
    pass


class Orbit(object):
    """Position and velocity of a body with respect to an attractor
    at a given time (epoch).

    """

    def __init__(self, state, epoch):
        """Constructor.

        Parameters
        ----------
        state : BaseState
            Position and velocity or orbital elements.
        epoch : ~astropy.time.Time
            Epoch of the orbit.

        """
        self._state = state  # type: BaseState
        self._epoch = epoch  # type: time.Time

    @property
    def state(self):
        """Position and velocity or orbital elements. """
        return self._state

    @property
    def epoch(self):
        """Epoch of the orbit. """
        return self._epoch

    @classmethod
    @u.quantity_input(r=u.m, v=u.m / u.s)
    def from_vectors(cls, attractor, r, v, epoch=J2000):
        """Return `Orbit` from position and velocity vectors.

        Parameters
        ----------
        attractor : Body
            Main attractor.
        r : ~astropy.units.Quantity
            Position vector wrt attractor center.
        v : ~astropy.units.Quantity
            Velocity vector.
        epoch : ~astropy.time.Time, optional
            Epoch, default to J2000.

        """
        assert np.any(r.value), "Position vector must be non zero"

        ss = rv.RVState(
            attractor, r, v)
        return cls(ss, epoch)

    @classmethod
    @u.quantity_input(a=u.m, ecc=u.one, inc=u.rad, raan=u.rad, argp=u.rad, nu=u.rad)
    def from_classical(cls, attractor, a, ecc, inc, raan, argp, nu, epoch=J2000):
        """Return `Orbit` from classical orbital elements.

        Parameters
        ----------
        attractor : Body
            Main attractor.
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
        epoch : ~astropy.time.Time, optional
            Epoch, default to J2000.

        """
        if ecc == 1.0 * u.one:
            raise ValueError("For parabolic orbits use Orbit.parabolic instead")

        if not 0 * u.deg <= inc <= 180 * u.deg:
            raise ValueError("Inclination must be between 0 and 180 degrees")

        if ecc > 1 and a > 0:
            raise ValueError("Hyperbolic orbits have negative semimajor axis")

        ss = classical.ClassicalState(
            attractor, a, ecc, inc, raan, argp, nu)
        return cls(ss, epoch)

    @classmethod
    @u.quantity_input(p=u.m, f=u.one, g=u.rad, h=u.rad, k=u.rad, L=u.rad)
    def from_equinoctial(cls, attractor, p, f, g, h, k, L, epoch=J2000):
        """Return `Orbit` from modified equinoctial elements.

        Parameters
        ----------
        attractor : Body
            Main attractor.
        p : ~astropy.units.Quantity
            Semilatus rectum.
        f : ~astropy.units.Quantity
            Second modified equinoctial element.
        g : ~astropy.units.Quantity
            Third modified equinoctial element.
        h : ~astropy.units.Quantity
            Fourth modified equinoctial element.
        k : ~astropy.units.Quantity
            Fifth modified equinoctial element.
        L : ~astropy.units.Quantity
            True longitude.
        epoch : ~astropy.time.Time, optional
            Epoch, default to J2000.

        """
        ss = equinoctial.ModifiedEquinoctialState(
            attractor, p, f, g, h, k, L)
        return cls(ss, epoch)

    @classmethod
    def from_body_ephem(cls, body, epoch=None):
        """Return osculating `Orbit` of a body at a given time.

        """
        if not epoch:
            epoch = time.Time.now().tdb
        elif epoch.scale != 'tdb':
            epoch = epoch.tdb
            warn("Input time was converted to scale='tdb' with value "
                 "{}. Use Time(..., scale='tdb') instead."
                 .format(epoch.tdb.value), TimeScaleWarning)

        r, v = get_body_barycentric_posvel(body.name, epoch)
        return cls.from_vectors(body.parent, r.xyz.to(u.km), v.xyz.to(u.km / u.day), epoch)

    @classmethod
    @u.quantity_input(alt=u.m, inc=u.rad, raan=u.rad, arglat=u.rad)
    def circular(cls, attractor, alt,
                 inc=0 * u.deg, raan=0 * u.deg, arglat=0 * u.deg, epoch=J2000):
        """Return circular `Orbit`.

        Parameters
        ----------
        attractor : Body
            Main attractor.
        alt : ~astropy.units.Quantity
            Altitude over surface.
        inc : ~astropy.units.Quantity, optional
            Inclination, default to 0 deg (equatorial orbit).
        raan : ~astropy.units.Quantity, optional
            Right ascension of the ascending node, default to 0 deg.
        arglat : ~astropy.units.Quantity, optional
            Argument of latitude, default to 0 deg.
        epoch: ~astropy.time.Time, optional
            Epoch, default to J2000.

        """
        a = attractor.R + alt
        ecc = 0 * u.one
        argp = 0 * u.deg

        return cls.from_classical(attractor, a, ecc, inc, raan, argp, arglat, epoch)

    @classmethod
    @u.quantity_input(p=u.m, inc=u.rad, raan=u.rad, argp=u.rad, nu=u.rad)
    def parabolic(cls, attractor, p, inc, raan, argp, nu, epoch=J2000):
        """Return parabolic `Orbit`.

        Parameters
        ----------
        attractor : Body
            Main attractor.
        p : ~astropy.units.Quantity
            Semilatus rectum or parameter.
        inc : ~astropy.units.Quantity, optional
            Inclination.
        raan : ~astropy.units.Quantity
            Right ascension of the ascending node.
        argp : ~astropy.units.Quantity
            Argument of the pericenter.
        nu : ~astropy.units.Quantity
            True anomaly.
        epoch: ~astropy.time.Time, optional
            Epoch, default to J2000.

        """
        k = attractor.k.to(u.km ** 3 / u.s ** 2)
        ecc = 1.0 * u.one
        r, v = classical.coe2rv(
            k.to(u.km ** 3 / u.s ** 2).value,
            p.to(u.km).value, ecc.value, inc.to(u.rad).value,
            raan.to(u.rad).value, argp.to(u.rad).value,
            nu.to(u.rad).value)

        ss = cls.from_vectors(attractor, r * u.km, v * u.km / u.s, epoch)
        return ss

    def __str__(self):
        if self.a > 1e7 * u.km:
            unit = u.au
        else:
            unit = u.km

        return ORBIT_FORMAT.format(
            r_p=self.r_p.to(unit).value, r_a=self.r_a.to(unit), inc=self.inc.to(u.deg),
            body=self.attractor
        )

    def __repr__(self):
        return self.__str__()

    def set_true_anomaly(self, target_nu):
        p, ecc, inc, raan, argp, _ = rv.rv2coe(self.attractor.k.to(u.km ** 3 / u.s ** 2).value,
                                               self.r.to(u.km).value,
                                               self.v.to(u.km / u.s).value)

        return self.from_classical(self.attractor, p / (1.0 - ecc ** 2) * u.km, 
                                   ecc * u.one, inc * u.rad, raan * u.rad, 
                                   argp * u.rad, target_nu * u.rad)

    def propagate(self, epoch_or_duration, function=propagate, rtol=1e-10):
        """Propagate this `Orbit` some `time` and return the result.

        Parameters
        ----------
        epoch_or_duration : Time, TimeDelta or equivalent
            Final epoch or time of flight.
        rtol : float, optional
            Relative tolerance for the propagation algorithm, default to 1e-10.

        """
        if isinstance(epoch_or_duration, time.Time) and not isinstance(epoch_or_duration, time.TimeDelta):
            time_of_flight = epoch_or_duration - self.epoch
        else:
            time_of_flight = time.TimeDelta(epoch_or_duration)

        return function(self, time_of_flight, rtol=rtol)

    def sample(self, values=None, function=propagate):
        """Samples an orbit to some specified time values.

        .. versionadded:: 0.8.0

        Parameters
        ----------
        values : Multiple options
            Number of interval points (default to 100),
            True anomaly values,
            Time values.

        Returns
        -------
        (Time, CartesianRepresentation)
            A tuple containing Time and Position vector in each
            given value.

        Notes
        -----
        When specifying a number of points, the initial and final
        position is present twice inside the result (first and
        last row). This is more useful for plotting.

        Examples
        --------
        >>> from astropy import units as u
        >>> from poliastro.examples import iss
        >>> iss.sample()
        >>> iss.sample(10)
        >>> iss.sample([0, 180] * u.deg)
        >>> iss.sample([0, 10, 20] * u.minute)
        >>> iss.sample([iss.epoch + iss.period / 2])

        """
        if values is None:
            return self.sample(100, function)

        elif isinstance(values, int):
            if self.ecc < 1:
                nu_values = np.linspace(0, 2 * np.pi, values) * u.rad
            else:
                # Select a sensible limiting value for non-closed orbits
                # This corresponds to r = 3p
                nu_limit = np.arccos(-(1 - 1 / 3.) / self.ecc)
                nu_values = np.linspace(-nu_limit, nu_limit, values)

            return self.sample(nu_values, function)

        elif hasattr(values, "unit") and values.unit in ('rad', 'deg'):
            values = self._generate_time_values(values)
        return (values, self._sample(values, function))

    def _sample(self, time_values, function=propagate):
        values = np.zeros((len(time_values), 3)) * self.r.unit
        for ii, epoch in enumerate(time_values):
            rr = self.propagate(epoch, function).r
            values[ii] = rr

        return CartesianRepresentation(values, xyz_axis=1)

    def _generate_time_values(self, nu_vals):
        M_vals = nu_to_M(nu_vals, self.ecc)
        time_values = self.epoch + (M_vals / self.n).decompose()
        return time_values

    def apply_maneuver(self, maneuver, intermediate=False):
        """Returns resulting `Orbit` after applying maneuver to self.

        Optionally return intermediate states (default to False).

        Parameters
        ----------
        maneuver : Maneuver
            Maneuver to apply.
        intermediate : bool, optional
            Return intermediate states, default to False.

        """
        orbit_new = self  # Initialize
        states = []
        attractor = self.attractor
        for delta_t, delta_v in maneuver:
            if not delta_t == 0 * u.s:
                orbit_new = orbit_new.propagate(delta_t)
            r, v = orbit_new.rv()
            vnew = v + delta_v
            orbit_new = self.from_vectors(attractor, r, vnew, orbit_new.epoch)
            states.append(orbit_new)
        if intermediate:
            res = states
        else:
            res = orbit_new
        return res

    # Delegated properties (syntactic sugar)
    def __getattr__(self, item):
        if hasattr(self.state, item):
            def delegated_(self_):
                return getattr(self_.state, item)

            # Use class docstring to properly translate properties, see
            # https://stackoverflow.com/a/38118315/554319
            delegated_.__doc__ = getattr(self.state.__class__, item).__doc__

            # Transform to a property
            delegated = property(delegated_)

        else:
            raise AttributeError("'{}' object has no attribute '{}'".format(self.__class__, item))

        # Bind the attribute
        setattr(self.__class__, item, delegated)

        # Return the newly bound attribute
        return getattr(self, item)
