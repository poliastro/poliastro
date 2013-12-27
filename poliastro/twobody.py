"""Two body problem.

"""

from astropy import units as u
from astropy import time

J2000 = time.Time("J2000", scale='utc')
ELEMENTS_UNITS = (u.m, u.dimensionless_unscaled, u.rad, u.rad, u.rad, u.rad)


class State(object):
    """Class to represent a position of a body wrt to an attractor.

    """
    def __init__(self, attractor, values, epoch=J2000):
        """Constructor.

        Parameters
        ----------
        attractor : Body
            Main attractor.
        values : tuple
            Values to define the state, either orbital elements
            or r, v vectors.
        epoch : Time, optional
            Epoch, default to J2000

        """
        self.attractor = attractor
        self.epoch = epoch
        if len(values) == 6:
            _check_elements_units(values)
            self.elements = values
        else:
            raise ValueError("Incorrect number of parameters")


def _check_elements_units(elements):
    """Check if orbital elements have consistent units.

    """
    for ii, unit in enumerate(ELEMENTS_UNITS):
        try:
            if elements[ii].si.unit != unit:
                raise u.UnitsError("Units must be consistent")
        except AttributeError:
            if ii != 1:
                raise ValueError("Elements must have units "
                                 "(use astropy.units)")
