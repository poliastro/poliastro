import numpy as np


class SpheroidLocation:
    """Class representing a ground station on an arbitrary ellipsoid."""

    def __init__(self, lon, lat, h, body):
        """

        Parameters
        ----------
        lon: ~astropy.units.quantity.Quantity
            geodetic longitude

        lat: ~astropy.units.quantity.Quantity
            geodetic latitude

        h: ~astropy.units.quantity.Quantity
            geodetic height

        body: ~poliastro.bodies.Body
            Planetary body the spheroid location lies on

        """
        self._lon = lon
        self._lat = lat
        self._h = h
        self._a = body.R
        self._b = body.R
        self._c = body.R_polar

    @property
    def cartesian_cords(self):
        """Convert to the Cartesian Coordinate system."""
        e2 = 1 - (self._c / self._a) ** 2
        N = self._a / np.sqrt(1 - e2 * np.sin(self._lon) ** 2)

        x = (N + self._h) * np.cos(self._lon) * np.cos(self._lat)
        y = (N + self._h) * np.cos(self._lon) * np.sin(self._lat)
        z = ((1 - e2) * N + self._h) * np.sin(self._lon)
        return x, y, z

    @property
    def f(self):
        """Get first flattening."""
        return 1 - self._c / self._a

    @property
    def N(self):
        """Normal vector of the ellipsoid at the given location."""
        x, y, z = self.cartesian_cords
        a, b, c = self._a.value, self._b.value, self._c.value
        N = np.array([2 * x.value / a ** 2, 2 * y.value / b ** 2, 2 * z.value / c ** 2])
        N /= np.linalg.norm(N)
        return N

    @property
    def tangential_vecs(self):
        """Returns orthonormal vectors tangential to the ellipsoid at the given location."""
        N = self.N
        u = np.array([1.0, 0, 0])
        u -= u.dot(N) * N
        u /= np.linalg.norm(u)
        v = np.cross(N, u)
        return u, v

    @property
    def radius_of_curvature(self):
        """Radius of curvature of the meridian at the latitude of the given location."""
        e2 = 1 - (self._c / self._a) ** 2
        rc = self._a * (1 - e2) / (1 - e2 * np.sin(self._lat) ** 2) ** 1.5
        return rc

    def distance(self, px, py, pz):
        """
        Calculates the distane from an arbitrary point to the given location (Cartesian coordinates).

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
        Determines whether an object located at a given point is visible from the given location.
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

    def cartesian_to_ellipsoidal(self, x, y, z):
        """
        Converts ellipsoidal coordinates to the Cartesian coordinate system for the given ellipsoid.
        Instead of the iterative formula, the function uses the approximation introduced in
        Bowring, B. R. (1976). TRANSFORMATION FROM SPATIAL TO GEOGRAPHICAL COORDINATES


        Parameters
        ----------
        x: ~astropy.units.quantity.Quantity
            x coordinate

        y: ~astropy.units.quantity.Quantity
            y coordinate

        z: ~astropy.units.quantity.Quantity
            z coordinate

        """
        a = self._a
        c = self._c
        e2 = 1 - (c / a) ** 2
        e2_ = e2 / (1 - e2)
        p = np.sqrt(x ** 2 + y ** 2)
        th = np.arctan(z * a / (p * c))
        lon = np.arctan(y / x)
        lat = np.arctan(
            (z + e2_ * c * np.sin(th) ** 3) / (p - e2 * a * np.cos(th) ** 3)
        )

        v = a / np.sqrt(1 - e2 * np.sin(lat) ** 2)
        h = p / np.cos(lat) - v
        h = z / np.sin(lat) - (1 - e2) * v

        return lat, lon, h
