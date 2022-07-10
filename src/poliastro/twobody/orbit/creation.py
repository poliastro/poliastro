from warnings import warn

import numpy as np
from astropy import units as u

from poliastro.constants import J2000
from poliastro.frames import Planes
from poliastro.frames.util import get_frame
from poliastro.twobody.elements import (
    get_eccentricity_critical_argp,
    get_eccentricity_critical_inc,
    get_inclination_critical_argp,
    heliosynchronous,
)
from poliastro.twobody.mean_elements import get_mean_elements
from poliastro.twobody.states import (
    ClassicalState,
    ModifiedEquinoctialState,
    RVState,
)
from poliastro.util import find_closest_value


class OrbitCreationMixin:
    """
    Mixin-class containing class-methods to create Orbit objects
    """

    def __init__(self, *_, **__):  # HACK stub to make mypy happy
        ...

    @classmethod
    @u.quantity_input(r=u.m, v=u.m / u.s)
    def from_vectors(
        cls, attractor, r, v, epoch=J2000, plane=Planes.EARTH_EQUATOR
    ):
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

        ss = RVState(attractor, (r, v), plane)
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
        attractor : Body
            Main attractor
        coord : ~astropy.coordinates.SkyCoord or ~astropy.coordinates.BaseCoordinateFrame
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
        inertial_frame_at_body_centre = get_frame(
            attractor, plane, coord.obstime
        )

        if not coord.is_equivalent_frame(inertial_frame_at_body_centre):
            coord_in_irf = coord.transform_to(inertial_frame_at_body_centre)
        else:
            coord_in_irf = coord

        pos = coord_in_irf.cartesian.xyz
        vel = coord_in_irf.cartesian.differentials["s"].d_xyz

        return cls.from_vectors(
            attractor, pos, vel, epoch=coord.obstime, plane=plane
        )

    @classmethod
    @u.quantity_input(
        a=u.m, ecc=u.one, inc=u.rad, raan=u.rad, argp=u.rad, nu=u.rad
    )
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
            raise ValueError(
                "For parabolic orbits use Orbit.parabolic instead"
            )

        if not 0 * u.deg <= inc <= 180 * u.deg:
            raise ValueError("Inclination must be between 0 and 180 degrees")

        if ecc > 1 and a > 0:
            raise ValueError("Hyperbolic orbits have negative semimajor axis")

        if not -np.pi * u.rad <= nu < np.pi * u.rad:
            warn("Wrapping true anomaly to -π <= nu < π", stacklevel=2)
            nu = (
                (nu + np.pi * u.rad) % (2 * np.pi * u.rad) - np.pi * u.rad
            ).to(nu.unit)

        ss = ClassicalState(
            attractor, (a * (1 - ecc**2), ecc, inc, raan, argp, nu), plane
        )
        return cls(ss, epoch)

    @classmethod
    @u.quantity_input(p=u.m, f=u.one, g=u.rad, h=u.rad, k=u.rad, L=u.rad)
    def from_equinoctial(
        cls,
        attractor,
        p,
        f,
        g,
        h,
        k,
        L,
        epoch=J2000,
        plane=Planes.EARTH_EQUATOR,
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
        ss = ModifiedEquinoctialState(attractor, (p, f, g, h, k, L), plane)
        return cls(ss, epoch)

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
        return cls.from_vectors(
            attractor, *ephem.rv(epoch), epoch, ephem.plane
        )

    @classmethod
    def from_sbdb(cls, name, **kwargs):
        """Return osculating `Orbit` by using `SBDB` from Astroquery.

        Parameters
        ----------
        name : str
            Name of the body to make the request.
        **kwargs
            Extra kwargs for astroquery.

        Returns
        -------
        ss: poliastro.twobody.orbit.Orbit
            Orbit corresponding to body_name

        Examples
        --------
        >>> from poliastro.twobody.orbit import Orbit
        >>> apophis_orbit = Orbit.from_sbdb('apophis')  # doctest: +REMOTE_DATA

        """
        from poliastro.io import orbit_from_sbdb

        return orbit_from_sbdb(name, **kwargs)

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
        epoch : ~astropy.time.Time, optional
            Epoch, default to J2000.
        plane : ~poliastro.frames.Planes
            Fundamental plane of the frame.

        """
        if alt < 0:
            raise ValueError("Altitude of an orbit cannot be negative.")
        a = attractor.R + alt
        ecc = 0 * u.one
        argp = 0 * u.deg

        return cls.from_classical(
            attractor=attractor,
            a=a,
            ecc=ecc,
            inc=inc,
            raan=raan,
            argp=argp,
            nu=arglat,
            epoch=epoch,
            plane=plane,
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
            raise ValueError(
                "The orbit for the given parameters doesn't exist"
            )

        return cls.from_classical(
            attractor=attractor,
            a=a_sync,
            ecc=ecc,
            inc=inc,
            raan=raan,
            argp=argp,
            nu=nu,
            epoch=epoch,
            plane=plane,
        )

    @classmethod
    def heliosynchronous(
        cls,
        attractor,
        a=None,
        ecc=None,
        inc=None,
        raan=0 * u.deg,
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
        attractor : ~poliastro.bodies.SolarSystemPlanet
            Attractor.
        a : ~astropy.units.Quantity
            Semi-major axis.
        ecc : ~astropy.units.Quantity
            Eccentricity.
        inc : ~astropy.units.Quantity
            Inclination.
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
        mean_elements = get_mean_elements(attractor)

        n_sunsync = (
            np.sqrt(mean_elements.attractor.k / abs(mean_elements.a**3))
            * u.one
        ).to(1 / u.s)

        try:
            a, ecc, inc = heliosynchronous(
                attractor.k, attractor.R, attractor.J2, n_sunsync, a, ecc, inc
            )
        except FloatingPointError:
            raise ValueError(
                "No SSO orbit with given parameters can be found."
            )

        ss = cls.from_classical(
            attractor=attractor,
            a=a,
            ecc=ecc,
            inc=inc,
            raan=raan,
            argp=argp,
            nu=nu,
            epoch=epoch.tdb,
            plane=plane,
        )

        return ss

    @classmethod
    @u.quantity_input(p=u.m, inc=u.rad, raan=u.rad, argp=u.rad, nu=u.rad)
    def parabolic(
        cls,
        attractor,
        p,
        inc,
        raan,
        argp,
        nu,
        epoch=J2000,
        plane=Planes.EARTH_EQUATOR,
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
        epoch : ~astropy.time.Time, optional
            Epoch, default to J2000.
        plane : ~poliastro.frames.Planes
            Fundamental plane of the frame.

        """
        ss = ClassicalState(
            attractor, (p, 1.0 * u.one, inc, raan, argp, nu), plane
        )
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
        epoch : ~astropy.time.Time, optional
            Epoch, default to J2000.
        plane : ~poliastro.frames.Planes
            Fundamental plane of the frame.

        """
        if attractor.J2 == 0.0 or attractor.J3 == 0.0:
            raise AttributeError(
                f"Attractor {attractor.name} has not spherical harmonics implemented"
            )

        critical_argps = [np.pi / 2, 3 * np.pi / 2] * u.rad
        critical_inclinations = [
            63.4349 * np.pi / 180,
            116.5651 * np.pi / 180,
        ] * u.rad
        try:
            if 1 <= np.abs(attractor.J2 / attractor.J3) <= 10:
                raise NotImplementedError(
                    f"This has not been implemented for {attractor.name}"
                )

            if alt < 0:
                raise ValueError("Altitude of an orbit cannot be negative")

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
                    attractor=attractor,
                    a=a,
                    ecc=ecc,
                    inc=inc,
                    raan=raan,
                    argp=argp,
                    nu=arglat,
                    epoch=epoch,
                    plane=plane,
                )

            inc = critical_inclinations[0] if inc is None else inc
            critical_inclination = find_closest_value(
                inc, critical_inclinations
            )
            if np.isclose(inc, critical_inclination, 1e-8, 1e-5 * u.rad):
                ecc = get_eccentricity_critical_inc(ecc)
                return cls.from_classical(
                    attractor=attractor,
                    a=a,
                    ecc=ecc,
                    inc=inc,
                    raan=raan,
                    argp=argp,
                    nu=arglat,
                    epoch=epoch,
                    plane=plane,
                )

            argp = critical_argps[0]
            inc = critical_inclinations[0]
            ecc = get_eccentricity_critical_argp(
                attractor.R, attractor.J2, attractor.J3, a, inc
            )

            return cls.from_classical(
                attractor=attractor,
                a=a,
                ecc=ecc,
                inc=inc,
                raan=raan,
                argp=argp,
                nu=arglat,
                epoch=epoch,
                plane=plane,
            )

        except AssertionError as exc:
            raise ValueError(
                f"The semimajor axis may not be smaller than the {attractor.name}'s radius"
            ) from exc
