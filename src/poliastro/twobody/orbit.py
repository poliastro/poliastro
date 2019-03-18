from warnings import warn

import numpy as np
from astropy import time, units as u
from astropy.coordinates import (
    GCRS,
    ICRS,
    Angle,
    CartesianDifferential,
    CartesianRepresentation,
    get_body_barycentric_posvel,
    solar_system_ephemeris,
)
from astroquery.jplhorizons import Horizons
from astroquery.jplsbdb import SBDB

from poliastro.bodies import Earth, Moon, Sun
from poliastro.constants import J2000
from poliastro.core.angles import nu_to_M as nu_to_M_fast
from poliastro.core.elements import rv2coe
from poliastro.frames import Planes, get_frame
from poliastro.plotting.core import OrbitPlotter2D, OrbitPlotter3D
from poliastro.twobody.angles import E_to_nu, M_to_nu, nu_to_M
from poliastro.twobody.propagation import mean_motion, propagate
from poliastro.util import hyp_nu_limit, norm

from ._states import BaseState, ClassicalState, ModifiedEquinoctialState, RVState

ORBIT_FORMAT = "{r_p:.0f} x {r_a:.0f} x {inc:.1f} ({frame}) orbit around {body} at epoch {epoch} ({scale})"
# String representation for orbits around bodies without predefined
# reference frame
ORBIT_NO_FRAME_FORMAT = (
    "{r_p:.0f} x {r_a:.0f} x {inc:.1f} orbit around {body} at epoch {epoch} ({scale})"
)


class TimeScaleWarning(UserWarning):
    pass


class OrbitSamplingWarning(UserWarning):
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
    def attractor(self):
        """Main attractor. """
        return self._state.attractor

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

    @property
    def r(self):
        """Position vector. """
        return self._state.to_vectors().r

    @property
    def v(self):
        """Velocity vector. """
        return self._state.to_vectors().v

    @property
    def a(self):
        """Semimajor axis. """
        return self.p / (1 - self.ecc ** 2)

    @property
    def p(self):
        """Semilatus rectum. """
        return self._state.to_classical().p

    @property
    def r_p(self):
        """Radius of pericenter. """
        return self.a * (1 - self.ecc)

    @property
    def r_a(self):
        """Radius of apocenter. """
        return self.a * (1 + self.ecc)

    @property
    def ecc(self):
        """Eccentricity. """
        return self._state.to_classical().ecc

    @property
    def inc(self):
        """Inclination. """
        return self._state.to_classical().inc

    @property
    def raan(self):
        """Right ascension of the ascending node. """
        return self._state.to_classical().raan

    @property
    def argp(self):
        """Argument of the perigee. """
        return self._state.to_classical().argp

    @property
    def nu(self):
        """True anomaly. """
        return self._state.to_classical().nu

    @property
    def f(self):
        """Second modified equinoctial element. """
        return self._state.to_equinoctial().f

    @property
    def g(self):
        """Third modified equinoctial element. """
        return self._state.to_equinoctial().g

    @property
    def h(self):
        """Fourth modified equinoctial element. """
        return self._state.to_equinoctial().h

    @property
    def k(self):
        """Fifth modified equinoctial element. """
        return self._state.to_equinoctial().k

    @property
    def L(self):
        """True longitude. """
        return self.raan + self.argp + self.nu

    @property
    def period(self):
        """Period of the orbit. """
        return 2 * np.pi * u.rad / self.n

    @property
    def n(self):
        """Mean motion. """
        return (np.sqrt(self.attractor.k / abs(self.a ** 3)) * u.rad).decompose()

    @property
    def energy(self):
        """Specific energy. """
        return self.v.dot(self.v) / 2 - self.attractor.k / np.sqrt(self.r.dot(self.r))

    @property
    def e_vec(self):
        """Eccentricity vector. """
        r, v = self.rv()
        k = self.attractor.k
        e_vec = ((v.dot(v) - k / (norm(r))) * r - r.dot(v) * v) / k
        return e_vec.decompose()

    @property
    def h_vec(self):
        """Specific angular momentum vector. """
        h_vec = np.cross(self.r.to(u.km).value, self.v.to(u.km / u.s)) * u.km ** 2 / u.s
        return h_vec

    @property
    def arglat(self):
        """Argument of latitude. """
        arglat = (self.argp + self.nu) % (360 * u.deg)
        return arglat

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

        ss = RVState(attractor, r, v)
        return cls(ss, epoch, plane)

    @classmethod
    def from_coords(cls, attractor, coord, plane=Planes.EARTH_EQUATOR):
        """Creates an `Orbit` from an attractor and astropy `SkyCoord`
        or `BaseCoordinateFrame` instance.

        This method accepts position
        and velocity in any reference frame unlike `Orbit.from_vector`
        which can accept inputs in only inertial reference frame centred
        at attractor. Also note that the frame information is lost after
        creation of the orbit and only the inertial reference frame at
        body centre will be used for all purposes.

        Parameters
        ----------
        attractor: Body
            Main attractor
        coord: astropy.coordinates.SkyCoord or BaseCoordinateFrame
            Position and velocity vectors in any reference frame. Note that coord must have
            a representation and its differential with respect to time.
        plane : ~poliastro.frames.Planes, optional
            Final orbit plane, default to Earth Equator.

        """
        if "s" not in coord.cartesian.differentials:
            raise ValueError(
                "Coordinate instance doesn't have a differential with respect to time"
            )
        if coord.size != 1:
            raise ValueError(
                "Coordinate instance must represents exactly 1 position, found: %d"
                % coord.size
            )

        # Reshape coordinate to 0 dimension if it is not already dimensionless.
        coord = coord.reshape(())

        # Get an inertial reference frame parallel to ICRS and centered at
        # attractor
        inertial_frame_at_body_centre = get_frame(attractor, plane, coord.obstime)

        if not coord.is_equivalent_frame(inertial_frame_at_body_centre):
            coord_in_irf = coord.transform_to(inertial_frame_at_body_centre)
        else:
            coord_in_irf = coord

        pos = coord_in_irf.cartesian.xyz
        vel = coord_in_irf.cartesian.differentials["s"].d_xyz

        return cls.from_vectors(attractor, pos, vel, epoch=coord.obstime, plane=plane)

    @classmethod
    @u.quantity_input(a=u.m, ecc=u.one, inc=u.rad, raan=u.rad, argp=u.rad, nu=u.rad)
    def from_classical(
        cls,
        attractor,
        a,
        ecc,
        inc,
        raan,
        argp,
        nu,
        epoch=J2000,
        plane=Planes.EARTH_EQUATOR,
    ):
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

        ss = ClassicalState(attractor, a * (1 - ecc ** 2), ecc, inc, raan, argp, nu)
        return cls(ss, epoch, plane)

    @classmethod
    @u.quantity_input(p=u.m, f=u.one, g=u.rad, h=u.rad, k=u.rad, L=u.rad)
    def from_equinoctial(
        cls, attractor, p, f, g, h, k, L, epoch=J2000, plane=Planes.EARTH_EQUATOR
    ):
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
        ss = ModifiedEquinoctialState(attractor, p, f, g, h, k, L)
        return cls(ss, epoch, plane)

    @classmethod
    def from_body_ephem(cls, body, epoch=None):

        """Return osculating `Orbit` of a body at a given time."""

        # TODO: https://github.com/poliastro/poliastro/issues/445

        if (
            body.name == "Pluto"
            and body.name.lower() not in solar_system_ephemeris.bodies
        ):
            raise KeyError(
                """These bodies can not be found in the "builtin" ephemeris. Tochange the ephemeri,. do
                >>> solar_system_ephemeris.set('jpl')"""
            )

        if not epoch:
            epoch = time.Time.now().tdb

        elif epoch.scale != "tdb":
            epoch = epoch.tdb
            warn(
                "Input time was converted to scale='tdb' with value "
                "{}. Use Time(..., scale='tdb') instead.".format(epoch.tdb.value),
                TimeScaleWarning,
            )

        r, v = get_body_barycentric_posvel(body.name, epoch)

        if body == Moon:
            # TODO: The attractor is in fact the Earth-Moon Barycenter
            icrs_cart = r.with_differentials(v.represent_as(CartesianDifferential))
            gcrs_cart = (
                ICRS(icrs_cart)
                .transform_to(GCRS(obstime=epoch))
                .represent_as(CartesianRepresentation)
            )
            ss = cls.from_vectors(
                Earth,
                gcrs_cart.xyz.to(u.km),
                gcrs_cart.differentials["s"].d_xyz.to(u.km / u.day),
                epoch,
            )

        else:
            # TODO: The attractor is not really the Sun, but the Solar System
            # Barycenter
            ss = cls.from_vectors(Sun, r.xyz.to(u.km), v.xyz.to(u.km / u.day), epoch)
            ss._frame = ICRS()  # Hack!

        return ss

    @classmethod
    def from_horizons(
        cls, name, epoch=None, plane=Planes.EARTH_EQUATOR, id_type="smallbody"
    ):
        """Return osculating `Orbit` of a body using JPLHorizons module of Astroquery.

        Parameters
        ----------
        name : string
            Name of the body to query for.
        epoch : ~astropy.time.Time, optional
            Epoch, default to None.
        plane : ~poliastro.frames.Planes
            Fundamental plane of the frame.
        id_type : string, optional
            Use "smallbody" for Asteroids and Comets, and "majorbody"
            for Planets and Satellites.

        """
        if not epoch:
            epoch = time.Time.now()
        if plane == Planes.EARTH_EQUATOR:
            refplane = "earth"
        elif plane == Planes.EARTH_ECLIPTIC:
            refplane = "ecliptic"

        obj = Horizons(id=name, epochs=epoch.jd, id_type=id_type).elements(
            refplane=refplane
        )
        a = obj["a"][0] * u.au
        ecc = obj["e"][0] * u.one
        inc = obj["incl"][0] * u.deg
        raan = obj["Omega"][0] * u.deg
        argp = obj["w"][0] * u.deg
        nu = obj["nu"][0] * u.deg
        ss = cls.from_classical(
            Sun, a, ecc, inc, raan, argp, nu, epoch=epoch.tdb, plane=plane
        )
        return ss

    @classmethod
    def from_sbdb(cls, name, **kargs):
        """Return osculating `Orbit` by using `SBDB` from Astroquery.

        Parameters
        ----------
        body_name: string
            Name of the body to make the request.

        Returns
        -------
        ss: poliastro.twobody.orbit.Orbit
            Orbit corresponding to body_name

        Examples
        --------
        >>> from poliastro.twobody.orbit import Orbit
        >>> apophis_orbit = Orbit.from_sbdb('apophis')
        """

        obj = SBDB.query(name, full_precision=True, **kargs)

        a = obj["orbit"]["elements"]["a"].to(u.AU) * u.AU
        ecc = float(obj["orbit"]["elements"]["e"]) * u.one
        inc = obj["orbit"]["elements"]["i"].to(u.deg) * u.deg
        raan = obj["orbit"]["elements"]["om"].to(u.deg) * u.deg
        argp = obj["orbit"]["elements"]["w"].to(u.deg) * u.deg

        # Since JPL provides Mean Anomaly (M) we need to make
        # the conversion to the true anomaly (\nu)
        nu = M_to_nu(obj["orbit"]["elements"]["ma"].to(u.deg) * u.deg, ecc)

        epoch = time.Time(obj["orbit"]["epoch"].to(u.d), format="jd")

        ss = cls.from_classical(
            Sun,
            a,
            ecc,
            inc,
            raan,
            argp,
            nu,
            epoch=epoch.tdb,
            plane=Planes.EARTH_ECLIPTIC,
        )

        return ss

    @classmethod
    @u.quantity_input(alt=u.m, inc=u.rad, raan=u.rad, arglat=u.rad)
    def circular(
        cls,
        attractor,
        alt,
        inc=0 * u.deg,
        raan=0 * u.deg,
        arglat=0 * u.deg,
        epoch=J2000,
        plane=Planes.EARTH_EQUATOR,
    ):
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

        return cls.from_classical(
            attractor, a, ecc, inc, raan, argp, arglat, epoch, plane
        )

    @classmethod
    @u.quantity_input(angular_velocity=u.rad / u.s, period=u.s, hill_radius=u.m)
    def geostationary(
        cls, attractor, angular_velocity=None, period=None, hill_radius=None
    ):
        """Return the geostationary orbit for the given attractor and its rotational speed.

        Parameters
        ----------
        attractor : Body
            Main attractor.
        angular_velocity : ~astropy.units.Quantity
            Rotational angular velocity of the attractor.
        period : ~astropy.units.Quantity
            Attractor's rotational period, ignored if angular_velocity is passed.
        hill_radius : ~astropy.units.Quantity
            Radius of Hill sphere of the attractor (optional). Hill sphere radius(in
            contrast with Laplace's SOI) is used here to validate the stability of the
            geostationary orbit, that is to make sure that the orbital radius required
            for the geostationary orbit is not outside of the gravitational sphere of
            influence of the attractor.
            Hill SOI of parent(if exists) of the attractor is ignored if hill_radius is not provided.
        """

        if angular_velocity is None and period is None:
            raise ValueError(
                "At least one among angular_velocity or period must be passed"
            )

        if angular_velocity is None:
            angular_velocity = 2 * np.pi / period

        # Find out geostationary radius using r = cube_root(GM/(angular
        # velocity)^2)
        with u.set_enabled_equivalencies(u.dimensionless_angles()):
            geo_radius = np.cbrt(attractor.k / np.square(angular_velocity.to(1 / u.s)))

        if hill_radius is not None and geo_radius > hill_radius:
            raise ValueError(
                "Geostationary orbit for the given parameters doesn't exist"
            )

        altitude = geo_radius - attractor.R
        return cls.circular(attractor, altitude)

    @classmethod
    @u.quantity_input(p=u.m, inc=u.rad, raan=u.rad, argp=u.rad, nu=u.rad)
    def parabolic(
        cls, attractor, p, inc, raan, argp, nu, epoch=J2000, plane=Planes.EARTH_EQUATOR
    ):
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
        ss = ClassicalState(attractor, p, 1.0 * u.one, inc, raan, argp, nu)
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
        >>> from astropy.coordinates import CartesianRepresentation, SphericalRepresentation
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
        # See the propagate function for reasoning about the usage of a
        # protected method
        coords = self.frame._replicate(cartesian, representation_type="cartesian")

        return coords.represent_as(representation)

    def to_icrs(self):
        """Creates a new Orbit object with its coordinates transformed to ICRS.

        Notice that, strictly speaking, the center of ICRS is the Solar System Barycenter
        and not the Sun, and therefore these orbits cannot be propagated in the context
        of the two body problem. Therefore, this function exists merely for practical
        purposes.

        .. versionadded:: 0.11.0

        """
        coords = self.frame.realize_frame(self.represent_as(CartesianRepresentation))
        coords.representation_type = CartesianRepresentation

        icrs_cart = coords.transform_to(ICRS).represent_as(CartesianRepresentation)

        # TODO: The attractor is in fact the Solar System Barycenter
        ss = self.from_vectors(
            Sun, r=icrs_cart.xyz, v=icrs_cart.differentials["s"].d_xyz, epoch=self.epoch
        )
        ss._frame = ICRS()  # Hack!
        return ss

    def rv(self):
        """Position and velocity vectors. """
        return self.r, self.v

    def classical(self):
        """Classical orbital elements. """
        return (
            self.a,
            self.ecc,
            self.inc.to(u.deg),
            self.raan.to(u.deg),
            self.argp.to(u.deg),
            self.nu.to(u.deg),
        )

    def pqw(self):
        """Perifocal frame (PQW) vectors. """
        if self.ecc < 1e-8:
            if abs(self.inc.to(u.rad).value) > 1e-8:
                node = np.cross([0, 0, 1], self.h_vec) / norm(self.h_vec)
                p_vec = node / norm(node)  # Circular inclined
            else:
                p_vec = [1, 0, 0] * u.one  # Circular equatorial
        else:
            p_vec = self.e_vec / self.ecc
        w_vec = self.h_vec / norm(self.h_vec)
        q_vec = np.cross(w_vec, p_vec) * u.one
        return p_vec, q_vec, w_vec

    def __str__(self):
        if self.a > 1e7 * u.km:
            unit = u.au
        else:
            unit = u.km

        try:
            return ORBIT_FORMAT.format(
                r_p=self.r_p.to(unit).value,
                r_a=self.r_a.to(unit),
                inc=self.inc.to(u.deg),
                frame=self.frame.__class__.__name__,
                body=self.attractor,
                epoch=self.epoch,
                scale=self.epoch.scale.upper(),
            )
        except NotImplementedError:
            return ORBIT_NO_FRAME_FORMAT.format(
                r_p=self.r_p.to(unit).value,
                r_a=self.r_a.to(unit),
                inc=self.inc.to(u.deg),
                body=self.attractor,
                epoch=self.epoch,
                scale=self.epoch.scale.upper(),
            )

    def __repr__(self):
        return self.__str__()

    def propagate(self, value, method=mean_motion, rtol=1e-10, **kwargs):
        """Propagates an orbit a specified time.

        If value is true anomaly, propagate orbit to this anomaly and return the result.
        Otherwise, if time is provided, propagate this `Orbit` some `time` and return the result.

        Parameters
        ----------
        value : ~astropy.units.Quantity, ~astropy.time.Time, ~astropy.time.TimeDelta
            Scalar time to propagate.
        rtol : float, optional
            Relative tolerance for the propagation algorithm, default to 1e-10.
        method : function, optional
            Method used for propagation
        **kwargs
            parameters used in perturbation models

        Returns
        -------
        Orbit
            New orbit after propagation.

        """
        if isinstance(value, time.Time) and not isinstance(value, time.TimeDelta):
            time_of_flight = value - self.epoch
        else:
            # Works for both Quantity and TimeDelta objects
            time_of_flight = time.TimeDelta(value)

        # TODO: Create a from_coordinates method
        coords = propagate(self, time_of_flight, method=method, rtol=rtol, **kwargs)[0]

        # Even after indexing the result of propagate,
        # the frame obstime might have an array of times
        # See the propagate function for reasoning about the usage of a
        # protected method
        coords = coords._replicate(
            coords.data, representation_type="cartesian", obstime=coords.obstime[0]
        )
        return self.from_coords(self.attractor, coords, plane=self.plane)

    @u.quantity_input(value=u.rad)
    def propagate_to_anomaly(self, value):
        """Propagates an orbit to a specific true anomaly.

        Parameters
        ----------
        value : ~astropy.units.Quantity

        Returns
        -------
        Orbit
            Resulting orbit after propagation.

        """
        # TODO: Avoid repeating logic with mean_motion?
        p, ecc, inc, raan, argp, _ = rv2coe(
            self.attractor.k.to(u.km ** 3 / u.s ** 2).value,
            self.r.to(u.km).value,
            self.v.to(u.km / u.s).value,
        )

        # Compute time of flight for correct epoch
        M = nu_to_M(self.nu, self.ecc)
        new_M = nu_to_M(value, self.ecc)
        time_of_flight = Angle(new_M - M).wrap_at(360 * u.deg) / self.n

        return self.from_classical(
            self.attractor,
            p / (1.0 - ecc ** 2) * u.km,
            ecc * u.one,
            inc * u.rad,
            raan * u.rad,
            argp * u.rad,
            value,
            epoch=self.epoch + time_of_flight,
            plane=self.plane,
        )

    def _sample_closed(self, values, min_anomaly, max_anomaly):
        limits = [
            min_anomaly.to(u.rad).value if min_anomaly is not None else 0,
            max_anomaly.to(u.rad).value if max_anomaly is not None else 2 * np.pi,
        ] * u.rad

        # First sample eccentric anomaly, then transform into true anomaly
        # to minimize error in the apocenter, see
        # http://www.dtic.mil/dtic/tr/fulltext/u2/a605040.pdf
        # Start from pericenter
        E_values = np.linspace(*limits, values)
        nu_values = E_to_nu(E_values, self.ecc)
        return nu_values

    def _sample_open(self, values, min_anomaly, max_anomaly):
        # Select a sensible limiting value for non-closed orbits
        # This corresponds to max(r = 3p, r = self.r)
        # We have to wrap nu in [-180, 180) to compare it with the output of
        # the arc cosine, which is in the range [0, 180)
        # Start from -nu_limit
        wrapped_nu = Angle(self.nu).wrap_at(180 * u.deg)
        nu_limit = max(hyp_nu_limit(self.ecc, 3.0), abs(wrapped_nu)).to(u.rad).value

        limits = [
            min_anomaly.to(u.rad).value if min_anomaly is not None else -nu_limit,
            max_anomaly.to(u.rad).value if max_anomaly is not None else nu_limit,
        ] * u.rad  # type: u.Quantity

        # Now we check that none of the provided values
        # is outside of the hyperbolic range
        nu_max = hyp_nu_limit(self.ecc) - 1e-3 * u.rad  # Arbitrary delta
        if not Angle(limits).is_within_bounds(-nu_max, nu_max):
            warn("anomaly outside range, clipping", OrbitSamplingWarning)
            limits = limits.clip(-nu_max, nu_max)

        nu_values = np.linspace(*limits, values)
        return nu_values

    def sample(
        self, values=100, *, min_anomaly=None, max_anomaly=None, method=mean_motion
    ):
        r"""Samples an orbit to some specified time values.

        .. versionadded:: 0.8.0

        Parameters
        ----------
        values : int
            Number of interval points (default to 100).
        min_anomaly, max_anomaly : ~astropy.units.Quantity, optional
            Anomaly limits to sample the orbit.
            For elliptic orbits the default will be :math:`E \in \left[0, 2 \pi \right]`,
            and for hyperbolic orbits it will be :math:`\nu \in \left[-\nu_c, \nu_c \right]`,
            where :math:`\nu_c` is either the current true anomaly
            or a value that corresponds to :math:`r = 3p`.
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
        >>> iss.sample()  # doctest: +ELLIPSIS
        <GCRS Coordinate ...>
        >>> iss.sample(10)  # doctest: +ELLIPSIS
        <GCRS Coordinate ...>

        """
        if self.ecc < 1:
            nu_values = self._sample_closed(values, min_anomaly, max_anomaly)
        else:
            nu_values = self._sample_open(values, min_anomaly, max_anomaly)

        time_values = self._generate_time_values(nu_values)
        return propagate(self, time.TimeDelta(time_values), method=method)

    def _generate_time_values(self, nu_vals):
        # Subtract current anomaly to start from the desired point
        ecc = self.ecc.value
        nu = self.nu.to(u.rad).value

        M_vals = [
            nu_to_M_fast(nu_val, ecc) - nu_to_M_fast(nu, ecc)
            for nu_val in nu_vals.to(u.rad).value
        ] * u.rad

        time_values = (M_vals / self.n).decompose()
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

    def plot(self, label=None, use_3d=False):
        """Plots the orbit as an interactive widget.

        Parameters
        ----------
        label : str, optional
            Label for the orbit, defaults to empty.
        use_3d : bool, optional
            Produce a 3D plot, default to False.

        """
        if use_3d:
            return OrbitPlotter3D().plot(self, label=label)
        else:
            return OrbitPlotter2D().plot(self, label=label)
