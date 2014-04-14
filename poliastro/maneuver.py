# coding: utf-8
"""Orbital maneuvers.

"""

import numpy as np
from numpy.linalg import norm

from astropy import units as u
u.one = u.dimensionless_unscaled  # astropy #1980

from poliastro.util import check_units


class Maneuver(object):
    """Class to represent a Maneuver.

    """
    def __init__(self, *vals):
        """Constructor.

        Parameters
        ----------
        vals : list
            List of pairs (delta_time, delta_velocity)

        """
        self.delta_times, self.delta_velocities = zip(*vals)
        u_times = check_units(self.delta_times, (u.s,))
        u_dvs = check_units(self.delta_velocities, (u.m / u.s,))
        if not u_times or not u_dvs:
            raise u.UnitsError("Units must be consistent")
        try:
            arr_dvs = [len(dv) == 3 for dv in self.delta_velocities]
            if not all(arr_dvs):
                raise TypeError
        except TypeError:
            raise ValueError("Delta-V must be three dimensions vectors")

    def get_total_time(self):
        """Returns total time of the maneuver.

        """
        total_time = sum(self.delta_times, 0 * u.s)
        return total_time

    def get_total_cost(self):
        """Returns otal cost of the maneuver.

        """
        dvs = [np.sqrt(dv.dot(dv)) for dv in self.delta_velocities]
        return sum(dvs, 0 * u.km / u.s)

    @classmethod
    def hohmann(cls, ss_i, r_f):
        """Compute a Hohmann transfer between two circular orbits.

        """
        # TODO: Check if circular?
        r_i = ss_i.a.to(u.km)
        v_i = ss_i.v.to(u.km / u.s)
        k = ss_i.attractor.k.to(u.km ** 3 / u.s ** 2)
        dv_a, dv_b, t_trans = hohmann(k, r_i, r_f)
        dv_vec_a = dv_a * v_i.value / norm(v_i.value)
        dv_vec_b = dv_b * v_i.value / norm(v_i.value)
        return cls((0 * u.s, dv_vec_a), (t_trans, dv_vec_b))


def hohmann(k, r_i, r_f):
    """Compute a Hohmann transfer between two circular orbits.

    Parameters
    ----------
    k : float
        Standard gravitational parameter.
    r_i, r_f : float
        Initial and final radiuses.

    Returns
    -------
    dv_a, dv_b : float
        Initial and final impulses.
    t_trans : float
        Time of transfer.

    """
    R = r_f / r_i
    v_i = np.sqrt(k / r_i)
    dv_a = v_i * (np.sqrt(2 * R / (1 + R)) - 1)
    dv_b = v_i * (1 - np.sqrt(2 / (1 + R))) / np.sqrt(R)
    # TODO: Refactor: half-period
    t_trans = np.pi * np.sqrt((r_i * (1 + R) / 2) ** 3 / k)
    return dv_a, dv_b, t_trans


def bielliptic(k, r_i, r_b, r_f):
    """Compute a bi-elliptic transfer between two circular orbits.

    Parameters
    ----------
    k : Quantity
        Standard gravitational parameter.
    r_i : Quantity
        Initial radius.
    r_b : Quantity
        Intermediate radius.
    r_f : Quantity
        Final radius.

    """
    R = r_f / r_i
    Rs = r_b / r_i
    v_i = np.sqrt(k / r_i)
    # TODO: Duplicate, see above
    dv_a = v_i * (np.sqrt(2 * Rs / (1 + Rs)) - 1)
    dv_b = v_i * np.sqrt(2 / Rs) * (np.sqrt(1 / (1 + Rs / R)) -
                                    np.sqrt(1 / (1 + Rs)))
    dv_c = v_i * (np.sqrt(2 * Rs / (R + Rs)) - 1) / np.sqrt(R)
    t_trans1 = np.pi * np.sqrt((r_i * (1 + Rs) / 2) ** 3 / k)
    t_trans2 = np.pi * np.sqrt((r_i * (R + Rs) / 2) ** 3 / k)
    return dv_a, dv_b, dv_c, t_trans1, t_trans2
