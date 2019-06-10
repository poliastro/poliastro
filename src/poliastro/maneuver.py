"""Orbital maneuvers.

"""
import numpy as np
from astropy import units as u

from poliastro.core.elements import coe_rotation_matrix, rv_pqw
from poliastro.core.util import cross
from poliastro.iod.izzo import lambert as lambert_izzo
from poliastro.util import norm


class Maneuver(object):
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
        impulses : list
            List of pairs (delta_time, delta_velocity)
        """

        self.impulses = args
        # HACK: Change API or validation code
        _dts, _dvs = zip(*args)
        self._dts, self._dvs = self._initialize(
            [(_dt * u.one).value for _dt in _dts] * (_dts[0] * u.one).unit,
            [(_dv * u.one).value for _dv in _dvs] * (_dvs[0] * u.one).unit,
        )
        try:
            if not all(len(dv) == 3 for dv in self._dvs):
                raise TypeError
        except TypeError:
            raise ValueError("Delta-V must be three dimensions vectors")

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
        dv: np.array
            Velocity components of the impulse.
        """

        return cls((0 * u.s, dv))

    @classmethod
    def hohmann(cls, orbit_i, r_f):
        r"""Compute a Hohmann transfer between two circular orbits.

        By defining the relationship between orbit radius:

        .. math::
            a_{trans} = \frac{r_{i} + r_{f}}{2}

        The Hohmann maneuver velocities can be expressed as:

        .. math::
            \begin{align}
                \Delta v_{a} &= \sqrt{\frac{2\mu}{r_{i}} - \frac{\mu}{a_{trans}}} - v_{i}\\
                \Delta v_{b} &= \sqrt{\frac{\mu}{r_{f}}} - \sqrt{\frac{2\mu}{r_{f}} - \frac{\mu}{a_{trans}}}
            \end{align}

        The time that takes to complete the maneuver can be computed as:

        .. math::
            \tau_{trans} = \pi \sqrt{\frac{(a_{trans})^{3}}{\mu}}

        Parameters
        ----------
        orbit_i: poliastro.twobody.orbit.Orbit
            Initial orbit
        r_f: astropy.unit.Quantity
            Final altitude of the orbit
        """

        if orbit_i.nu is not 0 * u.deg:
            orbit_i = orbit_i.propagate_to_anomaly(0 * u.deg)

        # Initial orbit data
        k = orbit_i.attractor.k
        r_i = orbit_i.r
        v_i = orbit_i.v
        h_i = norm(cross(r_i.to(u.m).value, v_i.to(u.m / u.s).value) * u.m ** 2 / u.s)
        p_i = h_i ** 2 / k.to(u.m ** 3 / u.s ** 2)

        # Hohmann is defined always from the PQW frame, since it is the
        # natural plane of the orbit
        r_i, v_i = rv_pqw(
            k.to(u.m ** 3 / u.s ** 2).value,
            p_i.to(u.m).value,
            orbit_i.ecc.value,
            orbit_i.nu.to(u.rad).value,
        )

        # Now, we apply Hohmman maneuver
        r_i = norm(r_i * u.m)
        v_i = norm(v_i * u.m / u.s)
        a_trans = (r_i + r_f) / 2

        # This is the modulus of the velocities
        dv_a = np.sqrt(2 * k / r_i - k / a_trans) - v_i
        dv_b = np.sqrt(k / r_f) - np.sqrt(2 * k / r_f - k / a_trans)

        # Write them in PQW frame
        dv_a = np.array([0, dv_a.to(u.m / u.s).value, 0])
        dv_b = np.array([0, -dv_b.to(u.m / u.s).value, 0])

        # Transform to IJK frame
        rot_matrix = coe_rotation_matrix(
            orbit_i.inc.to(u.rad).value,
            orbit_i.raan.to(u.rad).value,
            orbit_i.argp.to(u.rad).value,
        )

        dv_a = (rot_matrix @ dv_a) * u.m / u.s
        dv_b = (rot_matrix @ dv_b) * u.m / u.s

        t_trans = np.pi * np.sqrt(a_trans ** 3 / k)

        return cls((0 * u.s, dv_a), (t_trans, dv_b))

    @classmethod
    def bielliptic(cls, orbit_i, r_b, r_f):
        r"""Compute a bielliptic transfer between two circular orbits.

        The bielliptic maneuver employs two Hohmann transfers, therefore two
        intermediate orbits are established. We define the different radius
        relationships as follows:

        .. math::
            \begin{align}
                a_{trans1} &= \frac{r_{i} + r_{b}}{2}\\
                a_{trans2} &= \frac{r_{b} + r_{f}}{2}\\
            \end{align}

        The increments in the velocity are:

        .. math::
            \begin{align}
                \Delta v_{a} &= \sqrt{\frac{2\mu}{r_{i}} - \frac{\mu}{a_{trans1}}} - v_{i}\\
                \Delta v_{b} &= \sqrt{\frac{2\mu}{r_{b}} - \frac{\mu}{a_{trans2}}} - \sqrt{\frac{2\mu}{r_{b}} - \frac{\mu}{a_trans{1}}}\\
                \Delta v_{c} &= \sqrt{\frac{\mu}{r_{f}}} - \sqrt{\frac{2\mu}{r_{f}} - \frac{\mu}{a_{trans2}}}\\
            \end{align}

        The time of flight for this maneuver is the addition of the time needed for both transition orbits, following the same formula as
        Hohmann:

        .. math::
            \begin{align}
                \tau_{trans1} &= \pi \sqrt{\frac{a_{trans1}^{3}}{\mu}}\\
                \tau_{trans2} &= \pi \sqrt{\frac{a_{trans2}^{3}}{\mu}}\\
            \end{align}

        Parameters
        ----------
        orbit_i: poliastro.twobody.orbit.Orbit
            Initial orbit
        r_b: astropy.unit.Quantity
            Altitude of the intermediate orbit
        r_f: astropy.unit.Quantity
            Final altitude of the orbit
        """
        if orbit_i.nu is not 0 * u.deg:
            orbit_i = orbit_i.propagate_to_anomaly(0 * u.deg)

        # Initial orbit data
        k = orbit_i.attractor.k
        r_i = orbit_i.r
        v_i = orbit_i.v
        h_i = norm(cross(r_i.to(u.m).value, v_i.to(u.m / u.s).value) * u.m ** 2 / u.s)
        p_i = h_i ** 2 / k.to(u.m ** 3 / u.s ** 2)

        # Bielliptic is defined always from the PQW frame, since it is the
        # natural plane of the orbit
        r_i, v_i = rv_pqw(
            k.to(u.m ** 3 / u.s ** 2).value,
            p_i.to(u.m).value,
            orbit_i.ecc.value,
            orbit_i.nu.to(u.rad).value,
        )

        # Define the transfer radius
        r_i = norm(r_i * u.m)
        v_i = norm(v_i * u.m / u.s)
        a_trans1 = (r_i + r_b) / 2
        a_trans2 = (r_b + r_f) / 2

        # Compute impulses
        dv_a = np.sqrt(2 * k / r_i - k / a_trans1) - v_i
        dv_b = np.sqrt(2 * k / r_b - k / a_trans2) - np.sqrt(2 * k / r_b - k / a_trans1)
        dv_c = np.sqrt(k / r_f) - np.sqrt(2 * k / r_f - k / a_trans2)

        # Write impulses in PQW frame
        dv_a = np.array([0, dv_a.to(u.m / u.s).value, 0])
        dv_b = np.array([0, -dv_b.to(u.m / u.s).value, 0])
        dv_c = np.array([0, dv_c.to(u.m / u.s).value, 0])

        rot_matrix = coe_rotation_matrix(
            orbit_i.inc.to(u.rad).value,
            orbit_i.raan.to(u.rad).value,
            orbit_i.argp.to(u.rad).value,
        )

        # Transform to IJK frame
        dv_a = (rot_matrix @ dv_a) * u.m / u.s
        dv_b = (rot_matrix @ dv_b) * u.m / u.s
        dv_c = (rot_matrix @ dv_c) * u.m / u.s

        # Compute time for maneuver
        t_trans1 = np.pi * np.sqrt(a_trans1 ** 3 / k)
        t_trans2 = np.pi * np.sqrt(a_trans2 ** 3 / k)

        return cls((0 * u.s, dv_a), (t_trans1, dv_b), (t_trans2, dv_c))

    @classmethod
    def lambert(cls, orbit_i, orbit_f, method=lambert_izzo, short=True, **kwargs):
        """ Computes lambert maneuver between two different points.

        Parameters
        ----------
        orbit_i: ~poliastro.twobody.Orbit
            Initial orbit
        orbit_f: ~poliastro.twobody.Orbit
            Final orbit
        method: function
            Method for solving Lambert's problem
        short: keyword, boolean
            Selects between short and long solution
        """

        # Get initial algorithm conditions
        k = orbit_i.attractor.k
        r_i = orbit_i.r
        r_f = orbit_f.r

        # Time of flight is solved by substracting both orbit epochs
        tof = orbit_f.epoch - orbit_i.epoch

        # Compute all posible solutions to the Lambert transfer
        sols = list(method(k, r_i, r_f, tof, **kwargs))

        # Return short or long solution
        if short:
            dv_a, dv_b = sols[0]
        else:
            dv_a, dv_b = sols[-1]

        return cls((0 * u.s, dv_a - orbit_i.v), (tof.to(u.s), orbit_f.v - dv_b))

    def get_total_time(self):
        """Returns total time of the maneuver.

        """
        total_time = sum(self._dts, 0 * u.s)
        return total_time

    def get_total_cost(self):
        """Returns total cost of the maneuver.

        """
        dvs = [norm(dv) for dv in self._dvs]
        return sum(dvs, 0 * u.km / u.s)
