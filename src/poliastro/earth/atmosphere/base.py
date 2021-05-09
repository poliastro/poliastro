""" Holds different classes to model atmospheric models """

from poliastro.earth.atmosphere.util import h_to_z, z_to_h


class COESA:
    """Class for U.S Standard Atmosphere models."""

    def __init__(self, *tables):
        """Constructor for Atmosphere instances.

        Parameters
        ----------
        tables: *tables
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
        alt: ~astropy.units.Quantity
            Altitude to be checked.
        r0: ~astropy.units.Quantity
            Attractor radius.
        geometric: bool
            If `True`, assumes geometric altitude kind.

        Returns
        -------
        z: ~astropy.units.Quantity
            Geometric altitude.
        h: ~astropy.units.Quantity
            Geopotential altitude.

        """

        # Get geometric and geopotential altitudes
        if geometric is True:
            z = alt
            h = z_to_h(z, r0)
        else:
            h = alt
            z = h_to_z(h, r0)

        # Assert in range
        if z < self.zb_levels[0] or z > self.zb_levels[-1]:
            raise ValueError(
                "Geometric altitude must be in range [{}, {}]".format(
                    self.zb_levels[0], self.zb_levels[-1]
                )
            )

        return z, h

    def _get_index(self, x, x_levels):
        """Finds element in list and returns index.

        Parameters
        ----------
        x: ~astropy.units.Quantity
            Element to be searched.
        x_levels: list
            List for searching.

        Returns
        -------
        i: int
            Index for the value.

        """

        for i, value in enumerate(x_levels):
            if i < len(x_levels) and x > value:
                pass
            elif x == value:
                return i
            else:
                return i - 1
