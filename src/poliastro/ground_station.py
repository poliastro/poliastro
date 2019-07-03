import numpy as np

from astropy import units as u


class GroundStation(object):
    """"""

    def __init__(self, lon, lat, h, a, b):
        """

        Parameters
        ----------
        lon: ~astropy.units.core.Unit
            geodetic longitude

        lat: ~astropy.units.core.Unit
            geodetic latitude

        h: ~astropy.units.core.Unit
            geodetic height
        """
        self._lon = lon
        self._lat = lat
        self._h = h
        self._a = a
        self._b = b

    @property
    def cartesian_cords(self):
        """Convert to the Cartesian Coordinate system. """
        e2 = 1 - (self._b / self._a) ** 2
        N = self._a / np.sqrt(1 - e2 * np.sin(self._lon) ** 2)

        x = (N + self._h) * np.cos(self._lon) * np.cos(self._lat)
        y = (N + self._h) * np.cos(self._lon) * np.sin(self._lat)
        z = ((1 - e2) * N + self._h) * np.sin(self._lon)
        return x, y, z
