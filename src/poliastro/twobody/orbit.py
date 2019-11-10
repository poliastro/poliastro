from typing import List, Union
from warnings import warn

import numpy as np
from astropy import time, units as u
from astropy.coordinates import (
    GCRS,
    ICRS,
    Angle,
    CartesianDifferential,
    CartesianRepresentation,
    get_body_barycentric,
    get_body_barycentric_posvel,
)
from astroquery.jplhorizons import Horizons
from astroquery.jplsbdb import SBDB

from poliastro.bodies import Earth, Moon, Sun
from poliastro.constants import J2000
from poliastro.core.angles import nu_to_M as nu_to_M_fast
from poliastro.frames import Planes, get_frame
from poliastro.threebody.soi import laplace_radius
from poliastro.twobody.angles import E_to_nu, M_to_nu, nu_to_M, raan_from_ltan
from poliastro.twobody.propagation import mean_motion, propagate
from poliastro.util import (
    find_closest_value,
    get_eccentricity_critical_argp,
    get_eccentricity_critical_inc,
    get_inclination_critical_argp,
    hyp_nu_limit,
    norm,
)

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


class PatchedConicsWarning(UserWarning):
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
    def h_mag(self):
        """Specific angular momentum. """
        h_mag = norm(self.h_vec)
        return h_mag

    @property
    def arglat(self):
        """Argument of latitude. """
        arglat = (self.argp + self.nu) % (360 * u.deg)
        return arglat

    @property
    def M(self):
        """Mean anomaly. """
        return nu_to_M(self.nu, self.ecc)

    @property
    def t_p(self):
        """Elapsed time since latest perifocal passage. """
        t_p = self.period * self.M / (360 * u.deg)
        return t_p

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

        if r.ndim > 1 or v.ndim > 1:
            raise ValueError(
                f"Vectors must have dimension 1, got {r.ndim} and {v.ndim}"
            )

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
        """Return osculating `Orbit` of a body at a given time.

        """
        # TODO: https://github.com/poliastro/poliastro/issues/445
        if not epoch:
            epoch = time.Time.now().tdb
        elif epoch.scale != "tdb":
            epoch = epoch.tdb
            warn(
                "Input time was converted to scale='tdb' with value "
                f"{epoch.tdb.value}. Use Time(..., scale='tdb') instead.",
                TimeScaleWarning,
            )
        try:
            r, v = get_body_barycentric_posvel(body.name, epoch)
        except KeyError as exc:
            raise RuntimeError(
                """To compute the position and velocity of the Moon and Pluto use the JPL ephemeris:

>>> from astropy.coordinates import solar_system_ephemeris
>>> solar_system_ephemeris.set('jpl')
"""
            ) from exc
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

    def change_attractor(self, new_attractor, force=False):
        """Changes orbit attractor.

        Only changes from attractor to parent or the other way around are allowed.

        Parameters
        ----------
        new_attractor: poliastro.bodies.Body
            Desired new attractor.
        force: boolean
            If `True`, changes attractor even if physically has no-sense.

        Returns
        -------
        ss: poliastro.twobody.orbit.Orbit
            Orbit with new attractor

        """
        if self.attractor == new_attractor:
            return self
        elif self.attractor == new_attractor.parent:  # "Sun -> Earth"
            r_soi = laplace_radius(new_attractor)
            distance = norm(
                self.r - get_body_barycentric(new_attractor.name, self.epoch).xyz
            )
        elif self.attractor.parent == new_attractor:  # "Earth -> Sun"
            r_soi = laplace_radius(self.attractor)
            distance = norm(self.r)
        else:
            raise ValueError("Cannot change to unrelated attractor")

        if distance > r_soi and not force:
            raise ValueError(
                "Orbit is out of new attractor's SOI. If required, use 'force=True'."
            )
        elif self.ecc < 1.0 and not force:
            raise ValueError("Orbit will never leave the SOI of its current attractor")
        else:
            warn("Leaving the SOI of the current attractor", PatchedConicsWarning)

        new_frame = get_frame(new_attractor, self.plane, obstime=self.epoch)
        coords = self.frame.realize_frame(
            self.represent_as(CartesianRepresentation, CartesianDifferential)
        )
        ss = Orbit.from_coords(new_attractor, coords.transform_to(new_frame))

        return ss

    @classmethod
    def from_horizons(
        cls,
        name,
        attractor,
        epoch=None,
        plane=Planes.EARTH_EQUATOR,
        id_type="smallbody",
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

        bodies_dict = {
            "sun": 10,
            "mercury": 199,
            "venus": 299,
            "earth": 399,
            "mars": 499,
            "jupiter": 599,
            "saturn": 699,
            "uranus": 799,
            "neptune": 899,
        }
        location = "500@{}".format(bodies_dict[attractor.name.lower()])

        obj = Horizons(
            id=name, location=location, epochs=epoch.jd, id_type=id_type
        ).elements(refplane=refplane)
        a = obj["a"][0] * u.au
        ecc = obj["e"][0] * u.one
        inc = obj["incl"][0] * u.deg
        raan = obj["Omega"][0] * u.deg
        argp = obj["w"][0] * u.deg
        nu = obj["nu"][0] * u.deg
        ss = cls.from_classical(
            attractor, a, ecc, inc, raan, argp, nu, epoch=epoch.tdb, plane=plane
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
        >>> apophis_orbit = Orbit.from_sbdb('apophis')  # doctest: +REMOTE_DATA
        """

        obj = SBDB.query(name, full_precision=True, **kargs)

        if "count" in obj:
            # no error till now ---> more than one object has been found
            # contains all the name of the objects
            objects_name = obj["list"]["name"]
            objects_name_in_str = (
                ""  # used to store them in string form each in new line
            )
            for i in objects_name:
                objects_name_in_str += i + "\n"

            raise ValueError(
                str(obj["count"]) + " different objects found: \n" + objects_name_in_str
            )

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
    def heliosynchronous(
        cls,
        a=None,
        ecc=None,
        inc=None,
        ltan=10.0 * u.hourangle,
        argp=0 * u.deg,
        nu=0 * u.deg,
        epoch=J2000,
        body=Earth,
        plane=Planes.EARTH_EQUATOR,
    ):
        r""" Solves for a Sun-Synchronous orbit. These orbits make use of the J2
        perturbation to precess in order to be always towards Sun. At least
        two parameters of the set {a, ecc, inc} are needed in order to solve
        for these kind of orbits. Relationships among them are given by:

        .. math::
            \begin{align}
                a &= \left (\frac{-3R_{\bigoplus}J_{2}\sqrt{\mu}\cos(i)}{2\dot{\Omega}(1-e^2)^2}  \right ) ^ {\frac{2}{7}}\\
                e &= \sqrt{1 - \sqrt{\frac{-3R_{\bigoplus}J_{2}\sqrt{\mu}cos(i)}{2a^{\frac{7}{2}}\dot{\Omega}}}}\\
                i &= \arccos{\left ( \frac{-2a^{\frac{7}{2}}\dot{\Omega}(1-e^2)^2}{3R_{\bigoplus}J_{2}\sqrt{\mu}} \right )}\\
            \end{align}

        Parameters
        ----------
        a: ~astropy.units.Quantity
            Semi-major axis.
        ecc: ~astropy.units.Quantity
            Eccentricity.
        inc: ~astropy.units.Quantity
            Inclination.
        ltan: ~astropy.units.Quantity
            Local time of the ascending node which will be translated to the Right ascension of the ascending node.
        argp : ~astropy.units.Quantity
            Argument of the pericenter.
        nu : ~astropy.units.Quantity
            True anomaly.
        epoch : ~astropy.time.Time, optional
            Epoch, default to J2000.
        plane : ~poliastro.frames.Planes
            Fundamental plane of the frame.
        """

        # Constants for Sun-Synchronours Orbits (SSO)
        n_sunsync = 1.991063853e-7 * u.one / u.s
        R_SSO = body.R
        k_SSO = body.k
        J2_SSO = body.J2

        try:
            with np.errstate(invalid="raise"):
                if (a is None) and (ecc is None) and (inc is None):
                    # We check sufficient number of parameters
                    raise ValueError(
                        "At least two parameters of the set {a, ecc, inc} are required."
                    )
                elif a is None and (ecc is not None) and (inc is not None):
                    # Semi-major axis is the unknown variable
                    a = (
                        -3
                        * R_SSO ** 2
                        * J2_SSO
                        * np.sqrt(k_SSO)
                        / (2 * n_sunsync * (1 - ecc ** 2) ** 2)
                        * np.cos(inc)
                    ) ** (2 / 7)
                elif ecc is None and (a is not None) and (inc is not None):
                    # Eccentricity is the unknown variable
                    _ecc_0 = np.sqrt(
                        -3
                        * R_SSO ** 2
                        * J2_SSO
                        * np.sqrt(k_SSO)
                        * np.cos(inc.to(u.rad))
                        / (2 * a ** (7 / 2) * n_sunsync)
                    )
                    ecc = np.sqrt(1 - _ecc_0)
                elif inc is None and (ecc is not None) and (a is not None):
                    # Inclination is the unknown variable
                    inc = np.arccos(
                        -2
                        * a ** (7 / 2)
                        * n_sunsync
                        * (1 - ecc ** 2) ** 2
                        / (3 * R_SSO ** 2 * J2_SSO * np.sqrt(k_SSO))
                    )
        except FloatingPointError:
            raise ValueError("No SSO orbit with given parameters can be found.")
        raan = raan_from_ltan(epoch, ltan)
        ss = cls.from_classical(
            Earth, a, ecc, inc, raan, argp, nu, epoch=epoch.tdb, plane=plane
        )

        return ss

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

    @classmethod
    @u.quantity_input(
        alt=u.m, inc=u.rad, argp=u.rad, raan=u.rad, arglat=u.rad, ecc=u.one
    )
    def frozen(
        cls,
        attractor,
        alt,
        inc=None,
        argp=None,
        raan=0 * u.deg,
        arglat=0 * u.deg,
        ecc=None,
        epoch=J2000,
        plane=Planes.EARTH_EQUATOR,
    ):
        r"""Return frozen Orbit. If any of the given arguments results in an impossibility, some values will be overwritten

        To achieve frozen orbit these two equations have to be set to zero.

        .. math::
            \dfrac {d\overline {e}}{dt}=\dfrac {-3\overline {n}J_{3}R^{3}_{E}\sin \left( \overline {i}\right) }{2a^{3}\left( 1-\overline {e}^{2}\right) ^{2}}\left( 1-\dfrac {5}{4}\sin ^{2}\overline {i}\right) \cos \overline {w}

        .. math::
            \dfrac {d\overline {\omega }}{dt}=\dfrac {3\overline {n}J_{2}R^{2}_{E}}{a^{2}\left( 1-\overline {e}^{2}\right) ^{2}}\left( 1-\dfrac {5}{4}\sin ^{2}\overline {i}\right) \left[ 1+\dfrac {J_{3}R_{E}}{2J_{2}\overline {a}\left( 1-\overline {e}^{2}\right) }\left( \dfrac {\sin ^{2}\overline {i}-\overline {e}\cos ^{2}\overline {i}}{\sin \overline {i}}\right) \dfrac {\sin \overline {w}}{\overline {e}}\right]

        The first approach would be to nullify following term to zero:

        .. math::
            ( 1-\dfrac {5}{4}\sin ^{2})

        For which one obtains the so-called critical inclinations: i = 63.4349 or 116.5651 degrees. To escape the inclination requirement, the argument of periapsis can be set to w = 90 or 270 degrees to nullify the second equation. Then, one should nullify the right-hand side of the first equation, which yields an expression that correlates the inclination of the object and the eccentricity of the orbit:

        .. math::
            \overline {e}=-\dfrac {J_{3}R_{E}}{2J_{2}\overline {a}\left( 1-\overline {e}^{2}\right) }\left( \dfrac {\sin ^{2}\overline {i}-\overline {e}\cos ^{2} \overline {i}}{\sin \overline {i}}\right)

        Assuming that e is negligible compared to J2, it can be shown that:

        .. math::
            \overline {e}\approx -\dfrac {J_{3}R_{E}}{2J_{2}\overline {a}}\sin \overline {i}

        The implementation is divided in the following cases:

            1. When the user gives a negative altitude, the method will raise a ValueError
            2. When the attractor has not defined J2 or J3, the method will raise an AttributeError
            3. When the attractor has J2/J3 outside of range 1 to 10 , the method will raise an NotImplementedError. Special case for Venus.See "Extension of the critical inclination" by Xiaodong Liu, Hexi Baoyin, and Xingrui Ma
            4. If argp is not given or the given argp is a critical value:

                * if eccentricity is none and inclination is none, the inclination is set with a critical value and the eccentricity is obtained from the last formula mentioned
                * if only eccentricity is none, we calculate this value with the last formula mentioned
                * if only inclination is none ,we calculate this value with the formula for eccentricity with critical argp.

            5. If inc is not given or the given inc is critical:

                * if the argp and the eccentricity is given we keep these values to create the orbit
                * if the eccentricity is given we keep this value, if not, default to the eccentricity of the Moon's orbit around the Earth

            6. if it's not possible to create an orbit with the the argp and the inclination given, both of them are set to the critical values and the eccentricity is calculate with the last formula


        Parameters
        ----------
        attractor : Body
            Main attractor.
        alt : ~astropy.units.Quantity
            Altitude over surface.
        inc : ~astropy.units.Quantity, optional
            Inclination, default to critical value.
        argp : ~astropy.units.Quantity, optional
            Argument of the pericenter, default to critical value.
        raan : ~astropy.units.Quantity, optional
            Right ascension of the ascending node, default to 0 deg.
        arglat : ~astropy.units.Quantity, optional
            Argument of latitude, default to 0 deg.
        ecc : ~astropy.units.Quantity
            Eccentricity, default to the eccentricity of the Moon's orbit around the Earth
        epoch: ~astropy.time.Time, optional
            Epoch, default to J2000.
        plane : ~poliastro.frames.Planes
            Fundamental plane of the frame.
        """
        if attractor.J2 == 0.0 or attractor.J3 == 0.0:
            raise AttributeError(
                f"Attractor {attractor.name} has not spherical harmonics implemented"
            )

        critical_argps = [np.pi / 2, 3 * np.pi / 2] * u.rad
        critical_inclinations = [63.4349 * np.pi / 180, 116.5651 * np.pi / 180] * u.rad
        try:
            if 1 <= np.abs(attractor.J2 / attractor.J3) <= 10:
                raise NotImplementedError(
                    f"This has not been implemented for {attractor.name}"
                )

            assert alt > 0

            argp = critical_argps[0] if argp is None else argp
            a = attractor.R + alt

            critical_argp = find_closest_value(argp, critical_argps)
            if np.isclose(argp, critical_argp, 1e-8, 1e-5 * u.rad):
                if inc is None and ecc is None:
                    inc = critical_inclinations[0]
                    ecc = get_eccentricity_critical_argp(attractor, a, inc)
                elif ecc is None:
                    ecc = get_eccentricity_critical_argp(attractor, a, inc)
                else:
                    inc = get_inclination_critical_argp(attractor, a, ecc)
                return cls.from_classical(
                    attractor, a, ecc, inc, raan, argp, arglat, epoch, plane
                )

            inc = critical_inclinations[0] if inc is None else inc
            critical_inclination = find_closest_value(inc, critical_inclinations)
            if np.isclose(inc, critical_inclination, 1e-8, 1e-5 * u.rad):
                ecc = get_eccentricity_critical_inc(ecc)
                return cls.from_classical(
                    attractor, a, ecc, inc, raan, argp, arglat, epoch, plane
                )

            argp = critical_argps[0]
            inc = critical_inclinations[0]
            ecc = get_eccentricity_critical_argp(attractor, a, inc)

            return cls.from_classical(
                attractor, a, ecc, inc, raan, argp, arglat, epoch, plane
            )

        except AssertionError as exc:
            raise ValueError(
                f"The semimajor axis may not be smaller that {attractor.name}'s radius"
            ) from exc

    def represent_as(self, representation, differential_class=None):
        """Converts the orbit to a specific representation.

        .. versionadded:: 0.11.0

        Parameters
        ----------
        representation : ~astropy.coordinates.BaseRepresentation
            Representation object to use. It must be a class, not an instance.
        differential_class : ~astropy.coordinates.BaseDifferential, optional
            Class in which the differential should be represented, default to None.

        Examples
        --------
        >>> from poliastro.examples import iss
        >>> from astropy.coordinates import SphericalRepresentation
        >>> iss.represent_as(CartesianRepresentation)
        <CartesianRepresentation (x, y, z) in km
            (859.07256, -4137.20368, 5295.56871)>
        >>> iss.represent_as(CartesianRepresentation).xyz
        <Quantity [  859.07256, -4137.20368,  5295.56871] km>
        >>> iss.represent_as(CartesianRepresentation, CartesianDifferential).differentials['s']
        <CartesianDifferential (d_x, d_y, d_z) in km / s
            (7.37289205, 2.08223573, 0.43999979)>
        >>> iss.represent_as(CartesianRepresentation, CartesianDifferential).differentials['s'].d_xyz
        <Quantity [7.37289205, 2.08223573, 0.43999979] km / s>
        >>> iss.represent_as(SphericalRepresentation, CartesianDifferential)
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

        return cartesian.represent_as(representation, differential_class)

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

        cartesian = propagate(self, time_of_flight, method=method, rtol=rtol, **kwargs)
        new_epoch = self.epoch + time_of_flight

        # If the frame supports obstime, set the time values
        try:
            kwargs = {}
            if "obstime" in self.frame.frame_attributes:
                kwargs["obstime"] = new_epoch
            else:
                warn(
                    f"Frame {self.frame.__class__} does not support 'obstime', time values were not returned"
                )

            # Use of a protected method instead of frame.realize_frame
            # because the latter does not let the user choose the representation type
            # in one line despite its parameter names, see
            # https://github.com/astropy/astropy/issues/7784
            coords = self.frame._replicate(
                cartesian, representation_type="cartesian", **kwargs
            )

            return self.from_coords(self.attractor, coords, plane=self.plane)

        except NotImplementedError:
            return self.from_vectors(
                self.attractor,
                cartesian[0].xyz,
                cartesian[0].differentials["s"].d_xyz,
                new_epoch,
            )

    @u.quantity_input(value=u.rad)
    def time_to_anomaly(self, value):
        """ Returns time required to be in a specific true anomaly.

        Parameters
        ----------
        value : ~astropy.units.Quantity

        Returns
        -------
        tof: ~astropy.units.Quantity
            Time of flight required.
        """

        # Compute time of flight for correct epoch
        M = nu_to_M(self.nu, self.ecc)
        new_M = nu_to_M(value, self.ecc)
        tof = Angle(new_M - M).wrap_at(360 * u.deg) / self.n

        return tof

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

        # Compute time of flight for correct epoch
        time_of_flight = self.time_to_anomaly(value)

        return self.from_classical(
            self.attractor,
            self.a,
            self.ecc,
            self.inc,
            self.raan,
            self.argp,
            value,
            epoch=self.epoch + time_of_flight,
            plane=self.plane,
        )

    def _sample_closed(self, values, min_anomaly, max_anomaly):
        min_anomaly = (
            min_anomaly.to(u.rad).value
            if min_anomaly is not None
            else self.nu.to(u.rad).value
        )
        max_anomaly = (
            max_anomaly.to(u.rad).value
            if max_anomaly is not None
            else self.nu.to(u.rad).value + 2 * np.pi
        )
        limits = [min_anomaly, max_anomaly] * u.rad

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

    def sample(self, values=100, *, min_anomaly=None, max_anomaly=None):
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
        <CartesianRepresentation (x, y, z) in km ...
        >>> iss.sample(10)  # doctest: +ELLIPSIS
        <CartesianRepresentation (x, y, z) in km ...

        """
        if self.ecc < 1:
            nu_values = self._sample_closed(values, min_anomaly, max_anomaly)
        else:
            nu_values = self._sample_open(values, min_anomaly, max_anomaly)

        time_values = time.TimeDelta(self._generate_time_values(nu_values))
        cartesian = propagate(self, time_values)

        return cartesian

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
            res = states  # type: Union[Orbit, List[Orbit]]
        else:
            res = orbit_new
        return res

    def plot(self, label=None, use_3d=False, interactive=False):
        """Plots the orbit as an interactive widget.

        Parameters
        ----------
        label : str, optional
            Label for the orbit, defaults to empty.
        use_3d : bool, optional
            Produce a 3D plot, default to False.
        interactive : bool, optional
            Produce an interactive (rather than static) image of the orbit, default to False.
            This option requires Plotly properly installed and configured for your environment.
        """
        if not interactive and use_3d:
            raise ValueError(
                "The static plotter does not support 3D, use `interactive=True`"
            )
        elif not interactive:
            from poliastro.plotting.static import StaticOrbitPlotter

            return StaticOrbitPlotter().plot(self, label=label)
        elif use_3d:
            from poliastro.plotting.core import OrbitPlotter3D

            return OrbitPlotter3D().plot(self, label=label)
        else:
            from poliastro.plotting.core import OrbitPlotter2D

            return OrbitPlotter2D().plot(self, label=label)
