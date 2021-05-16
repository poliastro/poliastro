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

from poliastro.constants import J2000
from poliastro.frames import Planes
from poliastro.frames.util import get_frame
from poliastro.threebody.soi import laplace_radius
from poliastro.twobody.propagation import farnocchia, propagate

from ..core.elements import coe2rv_many
from ..core.propagation.farnocchia import delta_t_from_nu as delta_t_from_nu_fast
from ..util import find_closest_value, norm
from ..warnings import OrbitSamplingWarning, PatchedConicsWarning, TimeScaleWarning
from .angles import D_to_nu, E_to_nu, F_to_nu, M_to_D, M_to_E, M_to_F, raan_from_ltan
from .elements import (
    get_eccentricity_critical_argp,
    get_eccentricity_critical_inc,
    get_inclination_critical_argp,
    hyp_nu_limit,
)
from .mean_elements import get_mean_elements
from .sampling import sample_closed
from .states import BaseState, ClassicalState, ModifiedEquinoctialState, RVState

try:
    from functools import cached_property  # type: ignore
except ImportError:
    from cached_property import cached_property  # type: ignore


ORBIT_FORMAT = "{r_p:.0f} x {r_a:.0f} x {inc:.1f} ({frame}) orbit around {body} at epoch {epoch} ({scale})"
# String representation for orbits around bodies without predefined
# Reference frame
ORBIT_NO_FRAME_FORMAT = (
    "{r_p:.0f} x {r_a:.0f} x {inc:.1f} orbit around {body} at epoch {epoch} ({scale})"
)


class Orbit:
    """Position and velocity of a body with respect to an attractor
    at a given time (epoch).

    Regardless of how the Orbit is created, the implicit
    reference system is an inertial one. For the specific case
    of the Solar System, this can be assumed to be the
    International Celestial Reference System or ICRS.

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
        self._frame = None  # HACK: Only needed for Orbit.from_body_ephem

    @property
    def attractor(self):
        """Main attractor."""
        return self._state.attractor

    @property
    def epoch(self):
        """Epoch of the orbit."""
        return self._epoch

    @property
    def plane(self):
        """Fundamental plane of the frame."""
        return self._state.plane

    @cached_property
    def r(self):
        """Position vector."""
        return self._state.to_vectors().r

    @cached_property
    def v(self):
        """Velocity vector."""
        return self._state.to_vectors().v

    @cached_property
    def a(self):
        """Semimajor axis."""
        return self._state.to_classical().a

    @cached_property
    def p(self):
        """Semilatus rectum."""
        return self._state.to_classical().p

    @cached_property
    def r_p(self):
        """Radius of pericenter."""
        return self.a * (1 - self.ecc)

    @cached_property
    def r_a(self):
        """Radius of apocenter."""
        return self.a * (1 + self.ecc)

    @cached_property
    def ecc(self):
        """Eccentricity."""
        return self._state.to_classical().ecc

    @cached_property
    def inc(self):
        """Inclination."""
        return self._state.to_classical().inc

    @cached_property
    def raan(self):
        """Right ascension of the ascending node."""
        return self._state.to_classical().raan

    @cached_property
    def argp(self):
        """Argument of the perigee."""
        return self._state.to_classical().argp

    @property
    def nu(self):
        """True anomaly."""
        return self._state.to_classical().nu

    @cached_property
    def f(self):
        """Second modified equinoctial element."""
        return self._state.to_equinoctial().f

    @cached_property
    def g(self):
        """Third modified equinoctial element."""
        return self._state.to_equinoctial().g

    @cached_property
    def h(self):
        """Fourth modified equinoctial element."""
        return self._state.to_equinoctial().h

    @cached_property
    def k(self):
        """Fifth modified equinoctial element."""
        return self._state.to_equinoctial().k

    @cached_property
    def L(self):
        """True longitude."""
        return self.raan + self.argp + self.nu

    @cached_property
    def period(self):
        """Period of the orbit."""
        return self._state.period

    @cached_property
    def n(self):
        """Mean motion."""
        return self._state.n

    @cached_property
    def energy(self):
        """Specific energy."""
        return self.v.dot(self.v) / 2 - self.attractor.k / np.sqrt(self.r.dot(self.r))

    @cached_property
    def e_vec(self):
        """Eccentricity vector."""
        r, v = self.rv()
        k = self.attractor.k
        e_vec = ((v.dot(v) - k / (norm(r))) * r - r.dot(v) * v) / k
        return e_vec.decompose()

    @cached_property
    def h_vec(self):
        """Specific angular momentum vector."""
        h_vec = np.cross(self.r.to(u.km).value, self.v.to(u.km / u.s)) * u.km ** 2 / u.s
        return h_vec

    @cached_property
    def h_mag(self):
        """Specific angular momentum."""
        h_mag = norm(self.h_vec)
        return h_mag

    @cached_property
    def arglat(self):
        """Argument of latitude."""
        arglat = (self.argp + self.nu) % (360 * u.deg)
        return arglat

    @cached_property
    def t_p(self):
        """Elapsed time since latest perifocal passage."""
        t_p = (
            delta_t_from_nu_fast(
                self.nu.to_value(u.rad),
                self.ecc.value,
                self.attractor.k.to_value(u.km ** 3 / u.s ** 2),
                self.r_p.to_value(u.km),
            )
            * u.s
        )
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

        if r.ndim != 1 or v.ndim != 1:
            raise ValueError(
                f"Vectors must have dimension 1, got {r.ndim} and {v.ndim}"
            )

        ss = RVState(attractor, r, v, plane)
        return cls(ss, epoch)

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
        coord: ~astropy.coordinates.SkyCoord or ~astropy.coordinates.BaseCoordinateFrame
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
        for element in a, ecc, inc, raan, argp, nu, epoch:
            if not element.isscalar:
                raise ValueError(f"Elements must be scalar, got {element}")

        if ecc == 1.0 * u.one:
            raise ValueError("For parabolic orbits use Orbit.parabolic instead")

        if not 0 * u.deg <= inc <= 180 * u.deg:
            raise ValueError("Inclination must be between 0 and 180 degrees")

        if ecc > 1 and a > 0:
            raise ValueError("Hyperbolic orbits have negative semimajor axis")

        if not -np.pi * u.rad <= nu < np.pi * u.rad:
            warn("Wrapping true anomaly to -π <= nu < π", stacklevel=2)
            nu = ((nu + np.pi * u.rad) % (2 * np.pi * u.rad) - np.pi * u.rad).to(
                nu.unit
            )

        ss = ClassicalState(
            attractor, a * (1 - ecc ** 2), ecc, inc, raan, argp, nu, plane
        )
        return cls(ss, epoch)

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
        ss = ModifiedEquinoctialState(attractor, p, f, g, h, k, L, plane)
        return cls(ss, epoch)

    @classmethod
    def from_body_ephem(cls, body, epoch=None):
        """Return osculating `Orbit` of a body at a given time."""
        from poliastro.bodies import Earth, Moon, Sun

        warn(
            "Orbit.from_body_ephem is deprecated and will be removed in a future release, "
            "use Ephem.from_body instead",
            DeprecationWarning,
            stacklevel=2,
        )

        if not epoch:
            epoch = time.Time.now().tdb
        elif epoch.scale != "tdb":
            epoch = epoch.tdb
            warn(
                "Input time was converted to scale='tdb' with value "
                f"{epoch.tdb.value}. Use Time(..., scale='tdb') instead.",
                TimeScaleWarning,
                stacklevel=2,
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

    @classmethod
    def from_ephem(cls, attractor, ephem, epoch):
        """Create osculating orbit from ephemerides at a given epoch.

        This will assume that the `Ephem` coordinates
        are expressed with respect the given body.

        Parameters
        ----------
        ephem : ~poliastro.ephem.Ephem
            Ephemerides object to use.
        attractor : ~poliastro.bodies.Body
            Body to use as attractor.
        epoch : ~astropy.time.Time
            Epoch to retrieve the osculating orbit at.

        """
        return cls.from_vectors(attractor, *ephem.rv(epoch), epoch, ephem.plane)

    def get_frame(self):
        """Get equivalent reference frame of the orbit.

        .. versionadded:: 0.14.0

        """
        # HACK: Only needed for Orbit.from_body_ephem
        if self._frame is not None:
            return self._frame

        return get_frame(self.attractor, self.plane, self.epoch)

    def change_attractor(self, new_attractor, force=False):
        """Changes orbit attractor.

        Only changes from attractor to parent or the other way around are allowed.

        Parameters
        ----------
        new_attractor: poliastro.bodies.Body
            Desired new attractor.
        force: bool
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
            barycentric_position = get_body_barycentric(new_attractor.name, self.epoch)
            # Transforming new_attractor's frame into frame of attractor
            new_attractor_r = (
                ICRS(barycentric_position)
                .transform_to(self.get_frame())
                .represent_as(CartesianRepresentation)
                .xyz
            )
            distance = norm(self.r - new_attractor_r)
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
            warn(
                "Leaving the SOI of the current attractor",
                PatchedConicsWarning,
                stacklevel=2,
            )

        new_frame = get_frame(new_attractor, self.plane, obstime=self.epoch)
        coords = self.get_frame().realize_frame(
            self.represent_as(CartesianRepresentation, CartesianDifferential)
        )
        ss = Orbit.from_coords(new_attractor, coords.transform_to(new_frame))

        return ss

    def change_plane(self, plane):
        """Changes fundamental plane.

        Parameters
        ----------
        plane : ~poliastro.frames.Planes
            Fundamental plane of the frame.

        """
        if plane is self.plane:
            return self

        coords_orig = self.get_frame().realize_frame(
            self.represent_as(CartesianRepresentation, CartesianDifferential)
        )

        dest_frame = get_frame(self.attractor, plane, obstime=self.epoch)

        coords_dest = coords_orig.transform_to(dest_frame)
        coords_dest.representation_type = CartesianRepresentation

        return Orbit.from_coords(self.attractor, coords_dest, plane=plane)

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
        name : str
            Name of the body to query for.
        epoch : ~astropy.time.Time, optional
            Epoch, default to None.
        plane : ~poliastro.frames.Planes
            Fundamental plane of the frame.
        id_type : str, optional
            Use "smallbody" for Asteroids and Comets, and "majorbody"
            for Planets and Satellites.

        """
        warn(
            "Orbit.from_horizons is deprecated and will be removed in a future release, "
            "use Ephem.from_horizons instead",
            DeprecationWarning,
            stacklevel=2,
        )

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
        location = f"500@{bodies_dict[attractor.name.lower()]}"

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
    def from_sbdb(cls, name, **kwargs):
        """Return osculating `Orbit` by using `SBDB` from Astroquery.

        Parameters
        ----------
        name: str
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
        from poliastro.bodies import Sun

        obj = SBDB.query(name, full_precision=True, **kwargs)

        if "count" in obj:
            # No error till now ---> more than one object has been found
            # Contains all the name of the objects
            objects_name = obj["list"]["name"]
            objects_name_in_str = (
                ""  # Used to store them in string form each in new line
            )
            for i in objects_name:
                objects_name_in_str += i + "\n"

            raise ValueError(
                str(obj["count"]) + " different objects found: \n" + objects_name_in_str
            )

        if "object" not in obj.keys():
            raise ValueError("Object {} not found".format(name))

        a = obj["orbit"]["elements"]["a"]
        ecc = float(obj["orbit"]["elements"]["e"]) * u.one
        inc = obj["orbit"]["elements"]["i"]
        raan = obj["orbit"]["elements"]["om"]
        argp = obj["orbit"]["elements"]["w"]

        # Since JPL provides Mean Anomaly (M) we need to make
        # the conversion to the true anomaly (nu)
        M = obj["orbit"]["elements"]["ma"].to(u.rad)
        # NOTE: It is unclear how this conversion should happen,
        # see https://ssd-api.jpl.nasa.gov/doc/sbdb.html
        if ecc < 1:
            M = (M + np.pi * u.rad) % (2 * np.pi * u.rad) - np.pi * u.rad
            nu = E_to_nu(M_to_E(M, ecc), ecc)
        elif ecc == 1:
            nu = D_to_nu(M_to_D(M))
        else:
            nu = F_to_nu(M_to_F(M, ecc), ecc)

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
    def stationary(cls, attractor):
        """Return the stationary orbit for the given attractor and its rotational speed.

        Parameters
        ----------
        attractor : Body
            Main attractor.

        Returns
        -------
        Orbit
            New orbit.


        """
        return cls.synchronous(attractor)

    @classmethod
    @u.quantity_input(
        ecc=u.one,
        inc=u.deg,
        argp=u.deg,
        arglat=u.deg,
        raan=u.deg,
        period_mul=u.one,
    )
    def synchronous(
        cls,
        attractor,
        period_mul=1 * u.one,
        ecc=0 * u.one,
        inc=0 * u.deg,
        argp=0 * u.deg,
        arglat=0 * u.deg,
        raan=0 * u.deg,
        epoch=J2000,
        plane=Planes.EARTH_EQUATOR,
    ):
        r""" Returns an orbit where the orbital period equals the rotation rate
        of the orbited body.  The synchronous altitude for any central body can
        directly be obtained from Kepler's Third Law by setting the orbit period
        P\ :sub:`sync`, equal to the rotation period of the central body
        relative to the fixed stars D\ :sup:`*`. In order to obtain this, it's
        important to match orbital period with sidereal rotation period.

        Parameters
        ----------
        attractor : Body
            Main attractor.
        period_mul : ~astropy.units.Quantity
            Multiplier, default to 1 to indicate that the period of the body is
            equal to the sidereal rotational period of the body being orbited,
            0.5 a period equal to half the average rotational period of the body
            being orbited, indicates a semi-synchronous orbit.
        ecc : ~astropy.units.Quantity
            Eccentricity,default to 0 as a stationary orbit.
        inc : ~astropy.units.Quantity
            Inclination,default to 0 deg.
        raan : ~astropy.units.Quantity
            Right ascension of the ascending node,default to 0 deg.
        argp : ~astropy.units.Quantity
            Argument of the pericenter,default to 0 deg.
        arglat : ~astropy.units.Quantity, optional
            Argument of latitude, default to 0 deg.
        epoch : ~astropy.time.Time, optional
            Epoch, default to J2000.
        plane : ~poliastro.frames.Planes
            Fundamental plane of the frame.

        Returns
        -------
        Orbit
            New orbit.

        Raises
        ------
        ValueError
            If the pericenter is smaller than the attractor's radius.

        Notes
        -----

        Thus:

        .. math::

            P_{s y n c}=D^{*} \\

            a_{s y n c}=\left(\mu / 4 \pi^{2}\right)^{1 / 3}\left(D^{*}\right)^{2 / 3}\\

            H_{s y n c}=a_{s y n c} - R_{p l a n e t}\\

        """
        period_sync = attractor.rotational_period * period_mul
        a_sync = (attractor.k * (period_sync / (2 * np.pi)) ** 2) ** (1 / 3)
        nu = arglat - argp
        r_pericenter = (1 - ecc) * a_sync
        if r_pericenter < attractor.R:
            raise ValueError("The orbit for the given parameters doesn't exist")

        return cls.from_classical(
            attractor, a_sync, ecc, inc, raan, argp, nu, epoch, plane
        )

    @classmethod
    def heliosynchronous(
        cls,
        attractor,
        a=None,
        ecc=None,
        inc=None,
        ltan=10.0 * u.hourangle,
        argp=0 * u.deg,
        nu=0 * u.deg,
        epoch=J2000,
        plane=Planes.EARTH_EQUATOR,
    ):
        r"""Creates a Sun-Synchronous orbit.

        These orbits make use of the J2
        perturbation to precess in order to be always towards Sun. At least
        two parameters of the set {a, ecc, inc} are needed in order to solve
        for these kind of orbits.

        Relationships among them are given by:

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
        mean_elements = get_mean_elements(attractor)

        n_sunsync = (
            np.sqrt(mean_elements.attractor.k / abs(mean_elements.a ** 3)) * u.one
        ).decompose()
        R_SSO = attractor.R
        k_SSO = attractor.k
        J2_SSO = attractor.J2

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

        # Temporary fix: raan_from_ltan works only for Earth
        if attractor.name.lower() != "earth":
            raise NotImplementedError("Attractors other than Earth not supported yet")

        raan = raan_from_ltan(epoch, ltan)
        ss = cls.from_classical(
            attractor, a, ecc, inc, raan, argp, nu, epoch=epoch.tdb, plane=plane
        )

        return ss

    @classmethod
    @u.quantity_input(p=u.m, inc=u.rad, raan=u.rad, argp=u.rad, nu=u.rad)
    def parabolic(
        cls, attractor, p, inc, raan, argp, nu, epoch=J2000, plane=Planes.EARTH_EQUATOR
    ):
        """Return a parabolic `Orbit`.

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
        ss = ClassicalState(attractor, p, 1.0 * u.one, inc, raan, argp, nu, plane)
        return cls(ss, epoch)

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
        r"""Return a frozen Orbit. If any of the given arguments results in an impossibility, some values will be overwritten

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
                    ecc = get_eccentricity_critical_argp(
                        attractor.R, attractor.J2, attractor.J3, a, inc
                    )
                elif ecc is None:
                    ecc = get_eccentricity_critical_argp(
                        attractor.R, attractor.J2, attractor.J3, a, inc
                    )
                else:
                    inc = get_inclination_critical_argp(
                        attractor.R, attractor.J2, attractor.J3, a, ecc
                    )
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
            ecc = get_eccentricity_critical_argp(
                attractor.R, attractor.J2, attractor.J3, a, inc
            )

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

    def rv(self):
        """Position and velocity vectors."""
        return self.r, self.v

    def classical(self):
        """Classical orbital elements."""
        return (
            self.a,
            self.ecc,
            self.inc.to(u.deg),
            self.raan.to(u.deg),
            self.argp.to(u.deg),
            self.nu.to(u.deg),
        )

    def pqw(self):
        """Perifocal frame (PQW) vectors."""
        warn(
            "Orbit.pqw is deprecated and will be removed in a future release",
            DeprecationWarning,
            stacklevel=2,
        )

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
                frame=self.get_frame().__class__.__name__,
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

    def propagate(self, value, method=farnocchia, rtol=1e-10, **kwargs):
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

        return self.from_vectors(
            self.attractor,
            cartesian[0].xyz,
            cartesian[0].differentials["s"].d_xyz,
            new_epoch,
            plane=self.plane,
        )

    @u.quantity_input(value=u.rad)
    def time_to_anomaly(self, value):
        """Returns time required to be in a specific true anomaly.

        Parameters
        ----------
        value : ~astropy.units.Quantity

        Returns
        -------
        tof: ~astropy.units.Quantity
            Time of flight required.

        """
        # Silently wrap anomaly
        nu = (value + np.pi * u.rad) % (2 * np.pi * u.rad) - np.pi * u.rad

        delta_t = (
            delta_t_from_nu_fast(
                nu.to_value(u.rad),
                self.ecc.value,
                self.attractor.k.to_value(u.km ** 3 / u.s ** 2),
                self.r_p.to_value(u.km),
            )
            * u.s
        )
        tof = delta_t - self.t_p
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
        # Silently wrap anomaly
        nu = (value + np.pi * u.rad) % (2 * np.pi * u.rad) - np.pi * u.rad

        # Compute time of flight for correct epoch
        time_of_flight = self.time_to_anomaly(nu)

        if time_of_flight < 0:
            if self.ecc >= 1:
                raise ValueError("True anomaly {:.2f} not reachable".format(value))
            else:
                # For a closed orbit, instead of moving backwards
                # we need to do another revolution
                time_of_flight = self.period + time_of_flight

        return self.from_classical(
            self.attractor,
            self.a,
            self.ecc,
            self.inc,
            self.raan,
            self.argp,
            nu,
            epoch=self.epoch + time_of_flight,
            plane=self.plane,
        )

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
            warn("anomaly outside range, clipping", OrbitSamplingWarning, stacklevel=2)
            limits = limits.clip(-nu_max, nu_max)

        nu_values = np.linspace(*limits, values)  # type: ignore
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
        positions: ~astropy.coordinates.CartesianRepresentation
            Array of x, y, z positions.

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
            nu_values = sample_closed(
                min_anomaly if min_anomaly is not None else self.nu,
                self.ecc,
                max_anomaly,
                values,
            )
        else:
            nu_values = self._sample_open(values, min_anomaly, max_anomaly)

        n = nu_values.shape[0]
        rr, vv = coe2rv_many(
            np.full(n, self.attractor.k.to(u.m ** 3 / u.s ** 2).value),
            np.full(n, self.p.to(u.m).value),
            np.full(n, self.ecc.value),
            np.full(n, self.inc.to(u.rad).value),
            np.full(n, self.raan.to(u.rad).value),
            np.full(n, self.argp.to(u.rad).value),
            nu_values.to(u.rad).value,
        )

        # Add units
        rr = (rr << u.m).to(u.km)
        vv = (vv << (u.m / u.s)).to(u.km / u.s)

        cartesian = CartesianRepresentation(
            rr, differentials=CartesianDifferential(vv, xyz_axis=1), xyz_axis=1
        )

        return cartesian

    def _generate_time_values(self, nu_vals):
        # Subtract current anomaly to start from the desired point
        ecc = self.ecc.value
        k = self.attractor.k.to_value(u.km ** 3 / u.s ** 2)
        q = self.r_p.to_value(u.km)

        time_values = [
            delta_t_from_nu_fast(nu_val, ecc, k, q)
            for nu_val in nu_vals.to(u.rad).value
        ] * u.s - self.t_p
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
        """Plots the orbit.

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
