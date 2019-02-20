"""Orbital maneuvers.

"""
import numpy as np
from astropy import units as u

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
            R = \frac{r_{f}}{r_{i}}

        The Hohmann maneuver velocities can be expressed as:

        .. math::
            \begin{align}
                \Delta v_{a} &= v_{i}\left ( \sqrt{\frac{2R}{1+R}} - 1\right )\\
                \Delta v_{b} &= \frac{v_{i}}{\sqrt{R}}\left (\sqrt{\frac{2}{1+R}} - 1 \right )\\
            \end{align}

        The time that takes to compelte the maneuver can be computed as:

        .. math::
            \tau_{trans} = \pi \sqrt{\frac{(r{i} + r{f})^{3}}{8\mu}}

        Parameters
        ----------
        orbit_i: poliastro.twobody.orbit.Orbit
            Initial orbit
        r_f: astropy.unit.Quantity
            Final altitude of the orbit
        """

        if orbit_i.ecc == 0:

            r_i = orbit_i.a
            v_i = orbit_i.v
            k = orbit_i.attractor.k
            R = r_f / r_i
            dv_a = ((np.sqrt(2 * R / (1 + R)) - 1) * v_i).decompose()
            dv_b = (-(1 - np.sqrt(2 / (1 + R))) / np.sqrt(R) * v_i).decompose()
            t_trans = (np.pi * np.sqrt((r_i * (1 + R) / 2) ** 3 / k)).decompose()

            return cls((0 * u.s, dv_a), (t_trans, dv_b))

        else:
            raise ValueError(
                "Hohmann can only be applied to circular orbits (ecc == 0)"
            )

    @classmethod
    def bielliptic(cls, orbit_i, r_b, r_f):
        r"""Compute a bielliptic transfer between two circular orbits.

        The bielliptic maneuver employs two Hohmann transfers, therefore two
        intermediate orbits are stablished. We define the different radius
        relationships as follows:

        .. math::
            \begin{align}
                R &= \frac{r_{f}}{r_{i}}\\
                R_{s} &= \frac{r_{b}}{r_{i}}\\
            \end{align}

        The increments in the velocity are:

        .. math::
            \begin{align}
                \Delta v_{a} &= v_{i}\left ( \sqrt{\frac{2R{s}}{1+R{s}}} - 1\right )\\
                \Delta v_{b} &= v_{i}\sqrt{\frac{2}{R_{s}}}\left (\sqrt{\frac{1}{1+R_{s}}} -  \sqrt{\frac{1}{1 + \frac{R_{s}}{R}}} \right )\\
                \Delta v_{c} &= \frac{v_{i}}{\sqrt{R}}\left ( 1 - \sqrt{\frac{2R_{s}}{R+R_{s}}} \right )\\
            \end{align}

        The time of flight for this maneuver is the addition of the time needed for both transition orbits, following the same formula as
        Hohmann:

        .. math::
            \begin{align}
                \tau_{trans_{a}} &= \pi \sqrt{\frac{(r_{i}(1 + R_{s}))^{3}}{8\mu}}\\
                \tau_{trans_{b}} &= \pi \sqrt{\frac{(r_{i}(R + R_{s}))^{3}}{8\mu}}\\
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
        r_i = orbit_i.a
        v_i = orbit_i.v
        k = orbit_i.attractor.k
        R = r_f / r_i
        Rs = r_b / r_i
        dv_a = ((np.sqrt(2 * Rs / (1 + Rs)) - 1) * v_i).decompose()
        dv_b = (
            -np.sqrt(2 / Rs) * (np.sqrt(1 / (1 + Rs / R)) - np.sqrt(1 / (1 + Rs))) * v_i
        ).decompose()
        dv_c = (-(np.sqrt(2 * Rs / (R + Rs)) - 1) / np.sqrt(R) * v_i).decompose()
        t_trans1 = (np.pi * np.sqrt((r_i * (1 + Rs) / 2) ** 3 / k)).decompose()
        t_trans2 = (np.pi * np.sqrt((r_i * (R + Rs) / 2) ** 3 / k)).decompose()
        return cls((0 * u.s, dv_a), (t_trans1, dv_b), (t_trans2, dv_c))

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
