from warnings import warn

import numpy as np

from astropy import units as u
from astropy import time

from astropy.coordinates import (
    Angle,
    CartesianRepresentation, CartesianDifferential,
    get_body_barycentric_posvel,
    ICRS, GCRS
)

from poliastro.constants import J2000
from poliastro.twobody.angles import nu_to_M, E_to_nu
from poliastro.twobody.propagation import propagate, mean_motion
from poliastro.core.elements import rv2coe

from poliastro.twobody import rv
from poliastro.twobody import classical
from poliastro.twobody import equinoctial

from poliastro.bodies import Sun, Earth, Moon
from poliastro.frames import get_frame, Planes

from ._base import BaseState  # flake8: noqa


ORBIT_FORMAT = "{r_p:.0f} x {r_a:.0f} x {inc:.1f} ({frame}) orbit around {body}"


class TimeScaleWarning(UserWarning):
    pass


class Orbit(object):
    """Position and velocity of a body with respect to an attractor
    at a given time (epoch).

    Regardless of how the Orbit is created, the implicit
    reference system is an inertial one. For the specific case
    of the Solar System, this can be assumed to be the
    International Celestial Reference System or ICRS.

    """

    def __init__(self, state, epoch, plane):
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
        self._plane = plane
        self._frame = None

    @property
    def state(self):
        """Position and velocity or orbital elements. """
        return self._state

    @property
    def epoch(self):
        """Epoch of the orbit. """
        return self._epoch

    @property
    def plane(self):
        """Fundamental plane of the frame. """
        return self._plane

    @property
    def frame(self):
        """Reference frame of the orbit.

        .. versionadded:: 0.11.0

        """
        if self._frame is None:
            self._frame = get_frame(self.attractor, self._plane, self.epoch)

        return self._frame

    @classmethod
    @u.quantity_input(r=u.m, v=u.m / u.s)
    def from_vectors(cls, attractor, r, v, epoch=J2000, plane=Planes.EARTH_EQUATOR):
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
        plane : ~poliastro.frames.Planes
            Fundamental plane of the frame.

        """
        assert np.any(r.value), "Position vector must be non zero"

        ss = rv.RVState(
            attractor, r, v)
        return cls(ss, epoch, plane)

    @classmethod
    @u.quantity_input(a=u.m, ecc=u.one, inc=u.rad, raan=u.rad, argp=u.rad, nu=u.rad)
    def from_classical(cls, attractor, a, ecc, inc, raan, argp, nu, epoch=J2000, plane=Planes.EARTH_EQUATOR):
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
        plane : ~poliastro.frames.Planes
            Fundamental plane of the frame.

        """
        if ecc == 1.0 * u.one:
            raise ValueError("For parabolic orbits use Orbit.parabolic instead")

        if not 0 * u.deg <= inc <= 180 * u.deg:
            raise ValueError("Inclination must be between 0 and 180 degrees")

        if ecc > 1 and a > 0:
            raise ValueError("Hyperbolic orbits have negative semimajor axis")

        ss = classical.ClassicalState(
            attractor, a * (1 - ecc ** 2), ecc, inc, raan, argp, nu)
        return cls(ss, epoch, plane)

    @classmethod
    @u.quantity_input(p=u.m, f=u.one, g=u.rad, h=u.rad, k=u.rad, L=u.rad)
    def from_equinoctial(cls, attractor, p, f, g, h, k, L, epoch=J2000, plane=Planes.EARTH_EQUATOR):
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
        plane : ~poliastro.frames.Planes
            Fundamental plane of the frame.

        """
        ss = equinoctial.ModifiedEquinoctialState(
            attractor, p, f, g, h, k, L)
        return cls(ss, epoch, plane)

    @classmethod
    def from_body_ephem(cls, body, epoch=None):
        """Return osculating `Orbit` of a body at a given time.

        """
        # TODO: https://github.com/poliastro/poliastro/issues/445
        if not epoch:
            epoch = time.Time.now().tdb
        elif epoch.scale != 'tdb':
            epoch = epoch.tdb
            warn("Input time was converted to scale='tdb' with value "
                 "{}. Use Time(..., scale='tdb') instead."
                 .format(epoch.tdb.value), TimeScaleWarning)

        r, v = get_body_barycentric_posvel(body.name, epoch)

        if body == Moon:
            # The attractor is in fact the Earth-Moon Barycenter
            ss = cls.from_vectors(Earth, r.xyz.to(u.km), v.xyz.to(u.km / u.day), epoch)

        else:
            # The attractor is not really the Sun, but the Solar System Barycenter
            ss = cls.from_vectors(Sun, r.xyz.to(u.km), v.xyz.to(u.km / u.day), epoch)
            ss._frame = ICRS()  # Hack!

        return ss

    @classmethod
    @u.quantity_input(alt=u.m, inc=u.rad, raan=u.rad, arglat=u.rad)
    def circular(cls, attractor, alt,
                 inc=0 * u.deg, raan=0 * u.deg, arglat=0 * u.deg, epoch=J2000, plane=Planes.EARTH_EQUATOR):
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
        plane : ~poliastro.frames.Planes
            Fundamental plane of the frame.

        """
        a = attractor.R + alt
        ecc = 0 * u.one
        argp = 0 * u.deg

        return cls.from_classical(attractor, a, ecc, inc, raan, argp, arglat, epoch, plane)

    @classmethod
    @u.quantity_input(p=u.m, inc=u.rad, raan=u.rad, argp=u.rad, nu=u.rad)
    def parabolic(cls, attractor, p, inc, raan, argp, nu, epoch=J2000, plane=Planes.EARTH_EQUATOR):
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
        plane : ~poliastro.frames.Planes
            Fundamental plane of the frame.

        """
        ss = classical.ClassicalState(
            attractor, p, 1.0 * u.one, inc, raan, argp, nu)
        return cls(ss, epoch, plane)

    def represent_as(self, representation):
        """Converts the orbit to a specific representation.

        .. versionadded:: 0.11.0

        Parameters
        ----------
        representation : ~astropy.coordinates.BaseRepresentation
            Representation object to use. It must be a class, not an instance.

        Examples
        --------
        >>> from poliastro.examples import iss
        >>> from astropy.coordinates import CartesianRepresentation
        >>> iss.represent_as(CartesianRepresentation)
        <CartesianRepresentation (x, y, z) in km
            (859.07256, -4137.20368, 5295.56871)
         (has differentials w.r.t.: 's')>
        >>> iss.represent_as(CartesianRepresentation).xyz
        <Quantity [  859.07256, -4137.20368,  5295.56871] km>
        >>> iss.represent_as(CartesianRepresentation).differentials['s']
        <CartesianDifferential (d_x, d_y, d_z) in km / s
            (7.37289205, 2.08223573, 0.43999979)>
        >>> iss.represent_as(CartesianRepresentation).differentials['s'].d_xyz
        <Quantity [7.37289205, 2.08223573, 0.43999979] km / s>
        >>> iss.represent_as(SphericalRepresentation)
        <SphericalRepresentation (lon, lat, distance) in (rad, rad, km)
            (4.91712525, 0.89732339, 6774.76995296)
         (has differentials w.r.t.: 's')>

        """
        # As we do not know the differentials, we first convert to cartesian,
        # then let the frame represent_as do the rest
        # TODO: Perhaps this should be public API as well?
        cartesian = CartesianRepresentation(
            *self.r, differentials=CartesianDifferential(*self.v)
        )
        # See Orbit._sample for reasoning about the usage of a protected method
        coords = self.frame._replicate(cartesian, representation_type='cartesian')

        return coords.represent_as(representation)

    def __str__(self):
        if self.a > 1e7 * u.km:
            unit = u.au
        else:
            unit = u.km

        return ORBIT_FORMAT.format(
            r_p=self.r_p.to(unit).value, r_a=self.r_a.to(unit), inc=self.inc.to(u.deg),
            frame=self.frame.__class__.__name__,
            body=self.attractor,
        )

    def __repr__(self):
        return self.__str__()

    def propagate(self, value, method=mean_motion, rtol=1e-10, **kwargs):
        """Propagates an orbit.

        If value is true anomaly, propagate orbit to this anomaly and return the result.
        Otherwise, if time is provided, propagate this `Orbit` some `time` and return the result.

        Parameters
        ----------
        value : Multiple options
            True anomaly values or time values. If given an angle, it will always propagate forward.
        rtol : float, optional
            Relative tolerance for the propagation algorithm, default to 1e-10.
        method : function, optional
            Method used for propagation
        **kwargs
            parameters used in perturbation models

        """
        if hasattr(value, "unit") and value.unit in ('rad', 'deg'):
            p, ecc, inc, raan, argp, _ = rv2coe(self.attractor.k.to(u.km ** 3 / u.s ** 2).value,
                                                self.r.to(u.km).value,
                                                self.v.to(u.km / u.s).value)

            # Compute time of flight for correct epoch
            M = nu_to_M(self.nu, self.ecc)
            new_M = nu_to_M(value, self.ecc)
            time_of_flight = Angle(new_M - M).wrap_at(360 * u.deg) / self.n

            return self.from_classical(self.attractor, p / (1.0 - ecc ** 2) * u.km,
                                       ecc * u.one, inc * u.rad, raan * u.rad,
                                       argp * u.rad, value,
                                       epoch=self.epoch + time_of_flight, plane=self._plane)
        else:
            if isinstance(value, time.Time) and not isinstance(value, time.TimeDelta):
                time_of_flight = value - self.epoch
            else:
                time_of_flight = time.TimeDelta(value)

            return propagate(self, time_of_flight, method=method, rtol=rtol, **kwargs)

    def sample(self, values=None, method=mean_motion):
        """Samples an orbit to some specified time values.

        .. versionadded:: 0.8.0

        Parameters
        ----------
        values : Multiple options
            Number of interval points (default to 100),
            True anomaly values,
            Time values.

        method : function, optional
            Method used for propagation

        Returns
        -------
        positions: ~astropy.coordinates.BaseCoordinateFrame
            Array of x, y, z positions, with proper times as the frame attributes if supported.

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
            return self.sample(100, method)

        elif isinstance(values, int):
            if self.ecc < 1:
                # first sample eccentric anomaly, then transform into true anomaly
                # why sampling eccentric anomaly uniformly to minimize error in the apocenter, see
                # http://www.dtic.mil/dtic/tr/fulltext/u2/a605040.pdf
                # Start from pericenter
                E_values = np.linspace(0, 2 * np.pi, values) * u.rad
                nu_values = E_to_nu(E_values, self.ecc)
            else:
                # Select a sensible limiting value for non-closed orbits
                # This corresponds to max(r = 3p, r = self.r)
                # We have to wrap nu in [-180, 180) to compare it with the output of
                # the arc cosine, which is in the range [0, 180)
                # Start from -nu_limit
                wrapped_nu = self.nu if self.nu < 180 * u.deg else self.nu - 360 * u.deg
                nu_limit = max(np.arccos(-(1 - 1 / 3.) / self.ecc), abs(wrapped_nu))
                nu_values = np.linspace(-nu_limit, nu_limit, values)

            return self.sample(nu_values, method)

        elif hasattr(values, "unit") and values.unit in ('rad', 'deg'):
            values = self._generate_time_values(values)

        return self._sample(values, method)

    def _sample(self, time_values, method=mean_motion):
        positions = method(self, (time_values - self.epoch).to(u.s).value)

        data = CartesianRepresentation(positions[0] * u.km, xyz_axis=1)

        # If the frame supports obstime, set the time values
        kwargs = {}
        if 'obstime' in self.frame.frame_attributes:
            kwargs['obstime'] = time_values
        else:
            warn("Frame {} does not support 'obstime', time values were not returned".format(self.frame.__class__))

        # Use of a protected method instead of frame.realize_frame
        # because the latter does not let the user choose the representation type
        # in one line despite its parameter names, see
        # https://github.com/astropy/astropy/issues/7784
        return self.frame._replicate(data, representation_type='cartesian', **kwargs)

    def _generate_time_values(self, nu_vals):
        # Subtract current anomaly to start from the desired point
        M_vals = np.array([nu_to_M(nu_val, self.ecc).value -
                           nu_to_M(self.nu, self.ecc).value for nu_val in nu_vals]) * u.rad
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

    def __getstate__(self):
        return self.state, self.epoch

    def __setstate__(self, state):
        self._state = state[0]
        self._epoch = state[1]
