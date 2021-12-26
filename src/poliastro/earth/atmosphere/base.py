""" Holds different classes to model atmospheric models """

import astropy.units as u

from poliastro.core.earth_atmosphere.util import (
    _check_altitude as _check_altitude_fast,
    _get_index as _get_index_fast,
)


class COESA:
    """Class for U.S Standard Atmosphere models."""

    def __init__(self, *tables):
        """Constructor for Atmosphere instances.

        Parameters
        ----------
        *tables : *tables
            Different tables that define the atmosphere.

        """

        self.tables = tables

    @property
    def b_levels(self):
        return self.tables[0]

    @property
    def zb_levels(self):
        return self.tables[1]

    @property
    def hb_levels(self):
        return self.tables[2]

    @property
    def Tb_levels(self):
        return self.tables[3]

    @property
    def Lb_levels(self):
        return self.tables[4]

    @property
    def pb_levels(self):
        return self.tables[5]

    def _check_altitude(self, alt, r0, geometric=True):
        """Checks if altitude is inside valid range.

        Parameters
        ----------
        alt : ~astropy.units.Quantity
            Altitude to be checked.
        r0 : ~astropy.units.Quantity
            Attractor radius.
        geometric : bool
            If `True`, assumes geometric altitude kind.

        Returns
        -------
        z: ~astropy.units.Quantity
            Geometric altitude.
        h: ~astropy.units.Quantity
            Geopotential altitude.

        """
        alt = alt.to_value(u.km)
        r0 = r0.to_value(u.km)

        z, h = _check_altitude_fast(alt, r0, geometric)
        z, h = z * u.km, h * u.km

        # Assert in range
        if not self.zb_levels[0] <= z <= self.zb_levels[-1]:
            raise ValueError(
                f"Geometric altitude must be in range [{self.zb_levels[0]}, {self.zb_levels[-1]}]"
            )

        return z, h

    def _get_index(self, x, x_levels):
        """Finds element in list and returns index.

        Parameters
        ----------
        x : ~astropy.units.Quantity
            Element to be searched.
        x_levels : ~astropy.units.Quantity
            List for searching.

        Returns
        -------
        i: int
            Index for the value.

        """
        x = x.to_value(u.km)
        x_levels = (x_levels << u.km).value
        i = _get_index_fast(x, x_levels)
        return i
