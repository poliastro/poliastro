from typing import List, Union
from warnings import warn

import numpy as np
from astropy import time, units as u
from astropy.coordinates import (
    ICRS,
    CartesianDifferential,
    CartesianRepresentation,
    get_body_barycentric,
)

from poliastro.frames.util import get_frame
from poliastro.threebody.soi import laplace_radius
from poliastro.twobody.elements import (
    coe2rv_many,
    eccentricity_vector,
    energy,
    hyp_nu_limit,
    t_p,
)
from poliastro.twobody.propagation import farnocchia, propagate
from poliastro.twobody.sampling import sample_closed
from poliastro.twobody.states import BaseStateArray
from poliastro.util import norm, wrap_angle
from poliastro.warnings import OrbitSamplingWarning, PatchedConicsWarning

from .scalar import Orbit


ORBIT_FORMAT = "{r_p:.0f} x {r_a:.0f} x {inc:.1f} ({frame}) orbit around {body} at epoch {epoch} ({scale})"  # TODO import from orbit module?
# String representation for orbits around bodies without predefined
# Reference frame
ORBIT_NO_FRAME_FORMAT = (
    "{r_p:.0f} x {r_a:.0f} x {inc:.1f} orbit around {body} at epoch {epoch} ({scale})"
)  # TODO import from orbit module?


class OrbitArray:  # TODO creation mixin
    """Position and velocity of an array of body with respect to a common attractor
    at a given time (epoch).

    Regardless of how the Orbit is created, the implicit
    reference system is an inertial one. For the specific case
    of the Solar System, this can be assumed to be the
    International Celestial Reference System or ICRS.

    """

    def __init__(self, state):
        """Constructor.

        Parameters
        ----------
        state : BaseStateArray
            Position and velocity or orbital elements.

        """
        self._state = state  # type: BaseStateArray

    def __getitem__(self, idx):
        """Get item or slice from orbit array."""
        state = self._state[idx]
        cls = type(self) if isinstance(state, BaseStateArray) else Orbit
        return cls(self._state[idx])  # type: ignore

    def __setitem__(self, idx, value):
        """Set item or slice from array array."""
        raise NotImplementedError # TODO

    def copy(self):
        """Copy orbit array."""
        return type(self)(self._state.copy())

    def reshape(self, *args):
        """Reshape orbit array."""
        return type(self)(self._state.reshape(*args))

    @property
    def ndim(self):
        """Number of dimensions of array."""
        return self._state.ndim

    @property
    def shape(self):
        """Shape of array."""
        return self._state.shape

    @property
    def size(self):
        """Size of array."""
        return self._state.size

    @property
    def attractor(self):
        """Common main attractor."""
        return self._state.attractor

    @property
    def epoch(self):
        """Epoch array of the orbit array."""
        return self._state.epoch

    @property
    def plane(self):
        """Common fundamental plane of the frame."""
        return self._state.plane

    @property
    def r(self):
        """Position vector array."""
        return self._state.to_vectors().r

    @property
    def v(self):
        """Velocity vector array."""
        return self._state.to_vectors().v

    @property
    def a(self):
        """Semimajor axis array."""
        return self._state.to_classical().a

    @property
    def p(self):
        """Semilatus rectum array."""
        return self._state.to_classical().p

    @property
    def r_p(self):
        """Radius of pericenter array."""
        return self.a * (1 - self.ecc)

    @property
    def r_a(self):
        """Radius of apocenter array."""
        return self.a * (1 + self.ecc)

    @property
    def ecc(self):
        """Eccentricity array."""
        return self._state.to_classical().ecc

    @property
    def inc(self):
        """Inclination array."""
        return self._state.to_classical().inc

    @property
    def raan(self):
        """Right ascension of the ascending node array."""
        return self._state.to_classical().raan

    @property
    def argp(self):
        """Argument of the perigee array."""
        return self._state.to_classical().argp

    @property
    def nu(self):
        """True anomaly array."""
        return self._state.to_classical().nu

    @property
    def f(self):
        """Second modified equinoctial element array."""
        return self._state.to_equinoctial().f

    @property
    def g(self):
        """Third modified equinoctial element array."""
        return self._state.to_equinoctial().g

    @property
    def h(self):
        """Fourth modified equinoctial element array."""
        return self._state.to_equinoctial().h

    @property
    def k(self):
        """Fifth modified equinoctial element array."""
        return self._state.to_equinoctial().k

    @property
    def L(self):
        """True longitude array."""
        return self.raan + self.argp + self.nu

    @property
    def period(self):
        """Period of the orbit array."""
        return self._state.period

    @property
    def n(self):
        """Mean motion array."""
        return self._state.n

    @property
    def energy(self):
        """Specific energy array."""

        r = self.r  # on demand
        v = self.v  # on demand
        e = np.zeros(self.shape, dtype = r.dtype)
        k = self.attractor.k

        r_flat = r.reshape((-1, 3))
        v_flat = v.reshape((-1, 3))
        e_flat = e.reshape((-1,))

        for idx in range(self.size):  # TODO replace with vector implementation
            e_flat[idx] = energy(k, r_flat[idx, :], v_flat[idx, :])

        return e

    @property
    def e_vec(self):
        """Eccentricity vector array."""

        r = self.r  # on demand
        v = self.v  # on demand
        e = np.zeros((*self.shape, 3), dtype = r.dtype)
        k = self.attractor.k

        r_flat = r.reshape((-1, 3))
        v_flat = v.reshape((-1, 3))
        e_flat = e.reshape((-1, 3))

        for idx in range(self.size):  # TODO replace with vector implementation
            e_flat[idx, :] = eccentricity_vector(k, r_flat[idx, :], v_flat[idx, :])

        return e

    @property
    def h_vec(self):
        """Specific angular momentum vector array."""

        r = self.r  # on demand
        v = self.v  # on demand
        h_vec = np.zeros((*self.shape, 3), dtype = r.dtype)

        r_flat = r.to_value(u.km).reshape((-1, 3))
        v_flat = v.to(u.km / u.s).reshape((-1, 3))
        h_vec_flat = h_vec.reshape((-1, 3))

        for idx in range(self.size):  # TODO replace with vector implementation
            h_vec_flat[idx, :] = np.cross(r_flat[idx, :], v_flat[idx, :])

        return h_vec * u.km**2 / u.s

    @property
    def h_mag(self):
        """Specific angular momentum array."""

        h_vec = self.h_vec  # on demand
        h_mag = np.zeros(self.shape, dtype = h_vec.dtype)

        h_vec_flat = h_vec.reshape((-1, 3))
        h_mag_flat = h_mag.reshape((-1,))

        for idx in range(self.size):  # TODO replace with vector implementation
            h_mag_flat[idx] = norm(h_vec_flat[idx, :])

        return h_mag

    @property
    def arglat(self):
        """Argument of latitude array."""
        arglat = (self.argp + self.nu) % (360 * u.deg)
        return arglat

    @property
    def t_p(self):
        """Elapsed time since latest perifocal passage array."""

        nu = self.nu  # on demand
        ecc = self.ecc  # on demand
        k = self.attractor.k
        r_p = self.r_p  # on demand
        tp = np.zeros(self.shape, dtype = nu.dtype)

        nu_flat = nu.reshape((-1,))
        ecc_flat = ecc.reshape((-1,))
        r_p_flat = r_p.reshape((-1,))
        tp_flat = tp.reshape((-1,))

        for idx in range(self.size):  # TODO replace with vector implementation
            tp_flat[idx] = t_p(
                nu_flat[idx],
                ecc_flat[idx],
                k,
                r_p_flat[idx],
            )

        return tp

    def get_frame(self):
        """Get equivalent reference frame of the orbit.

        .. versionadded:: 0.14.0

        """
        # TODO returns object of subclass of `astropy.coordinates.BaseRADecFrame` - as dtype perhaps?
        # return get_frame(self.attractor, self.plane, self.epoch)
        raise NotImplementedError

    def change_attractor(self, new_attractor, force=False):
        """Changes orbit attractor.

        Only changes from attractor to parent or the other way around are allowed.

        Parameters
        ----------
        new_attractor : poliastro.bodies.Body
            Desired new attractor.
        force : bool
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
            if abs(self.inc.to_value(u.rad)) > 1e-8:
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
                r_p=self.r_p.to_value(unit),
                r_a=self.r_a.to(unit),
                inc=self.inc.to(u.deg),
                frame=self.get_frame().__class__.__name__,
                body=self.attractor,
                epoch=self.epoch,
                scale=self.epoch.scale.upper(),
            )
        except NotImplementedError:
            return ORBIT_NO_FRAME_FORMAT.format(
                r_p=self.r_p.to_value(unit),
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
        nu = wrap_angle(value)

        delta_t = t_p(
            nu,
            self.ecc,
            self.attractor.k,
            self.r_p,
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
                raise ValueError(f"True anomaly {value:.2f} not reachable")
            else:
                # For a closed orbit, instead of moving backwards
                # we need to do another revolution
                time_of_flight = self.period + time_of_flight

        return self.from_classical(
            attractor=self.attractor,
            a=self.a,
            ecc=self.ecc,
            inc=self.inc,
            raan=self.raan,
            argp=self.argp,
            nu=nu,
            epoch=self.epoch + time_of_flight,
            plane=self.plane,
        )

    def _sample_open(self, values, min_anomaly, max_anomaly):
        # Select a sensible limiting value for non-closed orbits
        # This corresponds to max(r = 3p, r = self.r)
        # We have to wrap nu in [-180, 180) to compare it with the output of
        # the arc cosine, which is in the range [0, 180)
        # Start from -nu_limit
        wrapped_nu = wrap_angle(self.nu)
        nu_limit = max(hyp_nu_limit(self.ecc, 3.0), abs(wrapped_nu)).to_value(u.rad)

        limits = [
            min_anomaly.to_value(u.rad) if min_anomaly is not None else -nu_limit,
            max_anomaly.to_value(u.rad) if max_anomaly is not None else nu_limit,
        ] * u.rad  # type: u.Quantity

        # Now we check that none of the provided values
        # is outside of the hyperbolic range
        nu_max = hyp_nu_limit(self.ecc) - 1e-3 * u.rad  # Arbitrary delta
        if not ((-nu_max <= limits).all() and (limits < nu_max).all()):
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
            np.tile(self.attractor.k, n),
            np.tile(self.p, n),
            np.tile(self.ecc, n),
            np.tile(self.inc, n),
            np.tile(self.raan, n),
            np.tile(self.argp, n),
            nu_values,
        )

        cartesian = CartesianRepresentation(
            rr, differentials=CartesianDifferential(vv, xyz_axis=1), xyz_axis=1
        )

        return cartesian

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
        plane = self.plane
        for delta_t, delta_v in maneuver:
            if not delta_t == 0 * u.s:
                orbit_new = orbit_new.propagate(delta_t)
            r, v = orbit_new.rv()
            vnew = v + delta_v
            orbit_new = self.from_vectors(
                attractor=attractor, r=r, v=vnew, epoch=orbit_new.epoch, plane=plane
            )
            if intermediate:  # Avoid keeping them in memory.
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
            from poliastro.plotting.interactive import OrbitPlotter3D

            return OrbitPlotter3D().plot(self, label=label)
        else:
            from poliastro.plotting.interactive import OrbitPlotter2D

            return OrbitPlotter2D().plot(self, label=label)
