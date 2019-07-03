import numpy as np


class GroundStation(object):
    """ Class representing a ground station on an arbitrary ellipsoid.
    """

    def __init__(self, lon, lat, h, a, b):
        """

        Parameters
        ----------
        lon: ~astropy.units.quantity.Quantity
            geodetic longitude

        lat: ~astropy.units.quantity.Quantity
            geodetic latitude

        h: ~astropy.units.quantity.Quantity
            geodetic height

        a: ~astropy.units.quantity.Quantity
            semi-major axis

        b: ~astropy.units.quantity.Quantity
            semi-minor axis
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

    @property
    def f(self):
        """Get first flattening. """
        return 1 - self._b / self._a

    @property
    def N(self):
        """Normal vector of the ellipsoid at the point of the ground station. """
        x, y, z = self.cartesian_cords
        N = np.array([2 * x.value, 2 * y.value, 2 * z.value])
        N = N / np.linalg.norm(N)
        return N

    @property
    def tangential_vecs(self):
        """Returns orthonormal vectors tangential to the ellipsoid at the point of the ground station. """
        N = self.N
        u = np.array([1.0, 0, 0])
        u -= u.dot(N) * N
        u /= np.linalg.norm(u)
        v = np.cross(N, u)
        return u, v

    def visible(self, px, py, pz):
        """
        Determines whether an object located at a given point is visible from the ground station.
        Returns true if true, false otherwise.

        Parameters
        ----------
        px: ~astropy.units.quantity.Quantity
            x-coordinate of the point
        py: ~astropy.units.quantity.Quantity
            y-coordinate of the point
        pz: ~astropy.units.quantity.Quantity
            z-coordinate of the point
        """
        c = [c.value for c in self.cartesian_cords]
        u = np.array([px.value, py.value, pz.value])

        N = self.N
        d = -N.dot(c)
        p = N.dot(u) + d

        return p >= 0
