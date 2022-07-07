"""Orbital maneuvers.

"""
from astropy import units as u

from poliastro.core.maneuver import (
    bielliptic as bielliptic_fast,
    correct_pericenter as correct_pericenter_fast,
    hohmann as hohmann_fast,
)
from poliastro.iod.izzo import lambert as lambert_izzo
from poliastro.util import norm


class Maneuver:
    r"""Class to represent a Maneuver.

    Each ``Maneuver`` consists on a list of impulses :math:`\Delta v_i`
    (changes in velocity) each one applied at a certain instant :math:`t_i`.
    You can access them directly indexing the ``Maneuver`` object itself.

    >>> man = Maneuver((0 * u.s, [1, 0, 0] * u.km / u.s),
    ... (10 * u.s, [1, 0, 0] * u.km / u.s))
    >>> man[0]
    (<Quantity 0. s>, <Quantity [1., 0., 0.] km / s>)
    >>> man.impulses[1]
    (<Quantity 10. s>, <Quantity [1., 0., 0.] km / s>)

    """

    def __init__(self, *args):
        r"""Constructor.

        Parameters
        ----------
        *args : list
            List of pairs (delta_time, delta_velocity)

        """

        self.impulses = args
        # HACK: Change API or validation code
        _dts, _dvs = zip(*args)

        try:
            self._dts, self._dvs = self._initialize(
                [(_dt * u.one).value for _dt in _dts] * (_dts[0] * u.one).unit,
                [(_dv * u.one).value for _dv in _dvs] * (_dvs[0] * u.one).unit,
            )
            if not all(len(dv) == 3 for dv in _dvs):
                raise TypeError
        except TypeError:
            raise ValueError("Delta-V must be three dimensions vectors")

    def __repr__(self):
        return f"Number of impulses: {len(self.impulses)}, Total cost: {self.get_total_cost():.6f}"

    @u.quantity_input(dts=u.s, dvs=u.m / u.s)
    def _initialize(self, dts, dvs):
        return dts, dvs

    def __getitem__(self, key):
        return self.impulses[key]

    @classmethod
    def impulse(cls, dv):
        """Single impulse at current time.

        Parameters
        ----------
        dv : numpy.ndarray
            Velocity components of the impulse.

        """

        return cls((0 * u.s, dv))

    @classmethod
    def hohmann(cls, orbit_i, r_f):
        r"""Compute a Hohmann transfer between two circular orbits.

        Parameters
        ----------
        orbit_i : poliastro.twobody.orbit.Orbit
            Initial orbit
        r_f : astropy.unit.Quantity
            Final orbital radius

        """

        # Propagate till periapsis
        if orbit_i.nu != 0 * u.deg:
            t_pericenter = orbit_i.time_to_anomaly(0 * u.deg)
            orbit_i = orbit_i.propagate_to_anomaly(0 * u.deg)
        else:
            t_pericenter = 0 * u.s

        k = orbit_i.attractor.k

        rv = orbit_i.rv()
        rv = (rv[0].to_value(u.m), rv[-1].to_value(u.m / u.s))

        k = k.to_value(u.m**3 / u.s**2)
        r_f = r_f.to_value(u.m)

        dv_a, dv_b, t_trans = hohmann_fast(k, rv, r_f)
        dv_a, dv_b, t_trans = dv_a * u.m / u.s, dv_b * u.m / u.s, t_trans * u.s

        return cls(
            (t_pericenter.decompose(), dv_a.decompose()),
            (t_trans.decompose(), dv_b.decompose()),
        )

    @classmethod
    def bielliptic(cls, orbit_i, r_b, r_f):
        r"""Compute a bielliptic transfer between two circular orbits.

        Parameters
        ----------
        orbit_i : poliastro.twobody.orbit.Orbit
            Initial orbit
        r_b : astropy.unit.Quantity
            Altitude of the intermediate orbit
        r_f : astropy.unit.Quantity
            Final orbital radius

        """

        # Propagate till periapsis
        if orbit_i.nu != 0 * u.deg:
            t_pericenter = orbit_i.time_to_anomaly(0 * u.deg)
            orbit_i = orbit_i.propagate_to_anomaly(0 * u.deg)
        else:
            t_pericenter = 0 * u.s

        k = orbit_i.attractor.k

        rv = orbit_i.rv()
        rv = (rv[0].to_value(u.m), rv[-1].to_value(u.m / u.s))

        k = k.to_value(u.m**3 / u.s**2)
        r_b = r_b.to_value(u.m)
        r_f = r_f.to_value(u.m)

        dv_a, dv_b, dv_c, t_trans1, t_trans2 = bielliptic_fast(
            k,
            r_b,
            r_f,
            rv,
        )
        dv_a, dv_b, dv_c, t_trans1, t_trans2 = (
            dv_a * u.m / u.s,
            dv_b * u.m / u.s,
            dv_c * u.m / u.s,
            t_trans1 * u.s,
            t_trans2 * u.s,
        )

        return cls(
            (t_pericenter.decompose(), dv_a.decompose()),
            (t_trans1.decompose(), dv_b.decompose()),
            (t_trans2.decompose(), dv_c.decompose()),
        )

    @classmethod
    def lambert(cls, orbit_i, orbit_f, method=lambert_izzo, **kwargs):
        """Computes Lambert maneuver between two different points.

        Parameters
        ----------
        orbit_i : ~poliastro.twobody.Orbit
            Initial orbit
        orbit_f : ~poliastro.twobody.Orbit
            Final orbit
        method : function
            Method for solving Lambert's problem
        **kwargs
            Extra kwargs for Lambert method.
        """

        # Get initial algorithm conditions
        k = orbit_i.attractor.k
        r_i = orbit_i.r
        r_f = orbit_f.r

        # Time of flight is solved by subtracting both orbit epochs
        tof = orbit_f.epoch - orbit_i.epoch
        if tof.to_value(u.s) < 0:
            raise ValueError(
                "Epoch of initial orbit greater than epoch of final orbit, "
                "causing a negative time of flight"
            )

        # Compute all possible solutions to the Lambert transfer
        dv_a, dv_b = method(k, r_i, r_f, tof, **kwargs)

        return cls(
            (0 * u.s, (dv_a - orbit_i.v).decompose()),
            (tof.to(u.s), (orbit_f.v - dv_b).decompose()),
        )

    def get_total_time(self):
        """Returns total time of the maneuver."""
        total_time = sum(self._dts, 0 * u.s)
        return total_time

    def get_total_cost(self):
        """Returns total cost of the maneuver."""
        dvs = [norm(dv) for dv in self._dvs]
        return sum(dvs, 0 * u.km / u.s)

    @classmethod
    @u.quantity_input(max_delta_r=u.km)
    def correct_pericenter(cls, orbit, max_delta_r):
        """Returns a Maneuver with the time before burning and the velocity vector in direction of the burn.

        Parameters
        ----------
        orbit : Orbit
            Position and velocity of a body with respect to an attractor
            at a given time (epoch).
        max_delta_r : ~astropy.units.Quantity
            Maximum satellite’s geocentric distance

        Returns
        -------
        maneuver: Maneuver
            Maneuver with the maximum time before we do an orbit-adjustment burn to restore the perigee to its
            nominal value and the velocity vector of the spacecraft to achieve the desired correction.

        Raises
        ------
        NotImplementedError
            - If the correction maneuver is not implemented for the attractor.
            - if the eccentricity is greater than 0.001.

        """
        J2 = orbit.attractor.J2.value
        if J2 == 0.0:
            raise NotImplementedError(
                f"The correction maneuver is not yet supported for {orbit.attractor}"
            )
        elif orbit.ecc > 0.001:
            raise NotImplementedError(
                f"The correction maneuver is not yet supported with {orbit.ecc},it should be less than or equal to 0.001"
            )

        R = orbit.attractor.R.to_value(u.km)
        k = orbit.attractor.k.to_value(u.km**3 / u.s**2)
        v = orbit.v.value
        a = orbit.a.value
        inc = orbit.inc.value
        ecc = orbit.ecc.value
        max_delta_r = max_delta_r.value

        delta_t, vf_ = correct_pericenter_fast(
            k, R, J2, max_delta_r, v, a, inc, ecc
        )
        delta_t = delta_t * u.s
        vf_ = vf_ * u.km / u.s

        return cls((delta_t, vf_))
