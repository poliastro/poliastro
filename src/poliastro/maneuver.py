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

    def __init__(self, *impulses):
        r"""Constructor.

        Parameters
        ----------
        impulses : list
            List of pairs (delta_time, delta_velocity)

        Notes
        -----
        TODO: Fix docstring, \*args convention

        """
        self.impulses = impulses
        # HACK: Change API or validation code
        _dts, _dvs = zip(*impulses)
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

        """
        return cls((0 * u.s, dv))

    @classmethod
    def hohmann(cls, orbit_i, r_f):
        """Compute a Hohmann transfer between two circular orbits.

        """
        # TODO: Check if circular?
        ss_i = orbit_i.state
        r_i = ss_i.a
        v_i = ss_i.v
        k = ss_i.attractor.k
        R = r_f / r_i
        dv_a = ((np.sqrt(2 * R / (1 + R)) - 1) * v_i).decompose()
        dv_b = (-(1 - np.sqrt(2 / (1 + R))) / np.sqrt(R) * v_i).decompose()
        t_trans = (np.pi * np.sqrt((r_i * (1 + R) / 2) ** 3 / k)).decompose()
        return cls((0 * u.s, dv_a), (t_trans, dv_b))

    @classmethod
    def bielliptic(cls, orbit_i, r_b, r_f):
        """Compute a bielliptic transfer between two circular orbits.

        """
        ss_i = orbit_i.state
        r_i = ss_i.a
        v_i = ss_i.v
        k = ss_i.attractor.k
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
