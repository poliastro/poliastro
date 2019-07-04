import numpy as np
from astropy import time, units as u


class GroundStation(object):
    """ Class representing a ground station on an arbitrary ellipsoid.
    """

    def __init__(self, lon, lat, h, epoch, angular_velocity, a, b, c=None):
        """

        Parameters
        ----------
        lon: ~astropy.units.quantity.Quantity
            geodetic longitude

        lat: ~astropy.units.quantity.Quantity
            geodetic latitude

        h: ~astropy.units.quantity.Quantity
            geodetic height

        epoch: ~astropy.units.quantity.Quantity
            starting epoch

        angular_velocity: ~astropy.units.quantity.Quantity
            angular velocity

        a: ~astropy.units.quantity.Quantity
            semi-major axis

        b: ~astropy.units.quantity.Quantity
            semi-minor axis

        c: ~astropy.units.quantity.Quantity
            Third axis defining the ellipsoid.
            If not set, it defaults to the value of a, making the ellipsoid
            a spheroid.
        """
        self._lon = lon
        self._lat = lat
        self._h = h
        self._epoch = epoch
        self._angular_velocity = angular_velocity
        self._a = a
        self._b = b
        self._c = c if c is not None else a

    @property
    def cartesian_cords(self):
        """Convert to the Cartesian Coordinate system.
        """
        e2 = 1 - (self._b / self._a) ** 2
        N = self._a / np.sqrt(1 - e2 * np.sin(self._lon) ** 2)

        x = (N + self._h) * np.cos(self._lon) * np.cos(self._lat)
        y = (N + self._h) * np.cos(self._lon) * np.sin(self._lat)
        z = ((1 - e2) * N + self._h) * np.sin(self._lon)
        return x, y, z

    @property
    def f(self):
        """Get first flattening.
        """
        return 1 - self._b / self._a

    @property
    def N(self):
        """Normal vector of the ellipsoid at the point of the ground station.
        """
        x, y, z = self.cartesian_cords
        a, b, c = self._a.value, self._b.value, self._c.value
        N = np.array([2 * x.value / a ** 2, 2 * y.value / b ** 2, 2 * z.value / c ** 2])
        N /= np.linalg.norm(N)
        return N

    @property
    def tangential_vecs(self):
        """Returns orthonormal vectors tangential to the ellipsoid at the point of the ground station.
        """
        N = self.N
        u = np.array([1.0, 0, 0])
        u -= u.dot(N) * N
        u /= np.linalg.norm(u)
        v = np.cross(N, u)
        return u, v

    def propagate(self, t):
        """
        Parameters
        ----------
        t: ~astropy.time.Time, ~astropy.time.TimeDelta
            Time to propagate, either an epoch
        """
        if isinstance(t, time.Time):
            delta_t = t - self._epoch
        elif isinstance(t, time.TimeDelta):
            delta_t = t
        else:
            raise ValueError(
                "t must be either of astropy.time.Time or astropy.time.TimeDelta type"
            )
        angle = delta_t * self._angular_velocity
        self._lat = (self._lat + angle) % (360 * u.deg)

    def distance(self, px, py, pz):
        """
        Calculates the distance from point to ground station.

        Parameters
        ----------
        px: ~astropy.units.quantity.Quantity
            x-coordinate of the point
        py: ~astropy.units.quantity.Quantity
            y-coordinate of the point
        pz: ~astropy.units.quantity.Quantity
            z-coordinate of the point
        """
        c = np.array([c.value for c in self.cartesian_cords])
        u = np.array([px.value, py.value, pz.value])
        d = np.linalg.norm(c - u)
        return d

    def is_visible(self, px, py, pz):
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
