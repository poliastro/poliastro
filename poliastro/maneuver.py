# coding: utf-8
"""Orbital maneuvers.

"""

import numpy as np

from astropy import units as u
u.one = u.dimensionless_unscaled  # astropy #1980

from poliastro.util import transform, check_units


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
        if not check_units(self.delta_times, (u.s,)) or not check_units(self.delta_velocities, (u.m / u.s,)):
            raise u.UnitsError("Units must be consistent")

    def total_time(self):
        """Total time of the maneuver.

        """
        # HACK: This won't work
        #total_time = sum(self.delta_times)
        total_time = 0.0 * u.s
        for dt in self.delta_times:
            total_time += dt
        return total_time

    @property
    def tof(self):
        """Total time of the maneuver.

        """
        return self.total_time()


def hohmann(k, r_i, r_f):
    """Compute a Hohmann transfer between two circular orbits.

    Parameters
    ----------
    k : Quantity
        Standard gravitational parameter.
    r_i, r_f : Quantity
        Initial and final radiuses.

    Returns
    -------
    dv_a, dv_b : Quantity
        Initial and final impulses.
    t_trans : Quantity
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
