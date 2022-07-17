import astropy.units as u
import numpy as np

from poliastro.core.spheroid_location import (
    N as N_fast,
    cartesian_cords as cartesian_cords_fast,
    cartesian_to_ellipsoidal as cartesian_to_ellipsoidal_fast,
    distance as distance_fast,
    f as f_fast,
    is_visible as is_visible_fast,
    radius_of_curvature as radius_of_curvature_fast,
    tangential_vecs as tangential_vecs_fast,
)


class SpheroidLocation:
    """Class representing a ground station on an oblate ellipsoid."""

    def __init__(self, lon, lat, h, body):
        """
        Parameters
        ----------
        lon : ~astropy.units.quantity.Quantity
            Geodetic longitude
        lat : ~astropy.units.quantity.Quantity
            Geodetic latitude
        h : ~astropy.units.quantity.Quantity
            Geodetic height
        body : ~poliastro.bodies.Body
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
        _a, _c, _lon, _lat, _h = (
            self._a.to_value(u.m),
            self._c.to_value(u.m),
            self._lon.to_value(u.rad),
            self._lat.to_value(u.rad),
            self._h.to_value(u.m),
        )
        return (
            cartesian_cords_fast(
                _a,
                _c,
                _lon,
                _lat,
                _h,
            )
            * u.m
        )

    @property
    def f(self):
        """Get first flattening."""
        _a, _c = self._a.to_value(u.m), self._c.to_value(u.m)
        return f_fast(_a, _c)

    @property
    def N(self):
        """Normal vector of the ellipsoid at the given location."""
        a, b, c = (
            self._a.to_value(u.m),
            self._b.to_value(u.m),
            self._c.to_value(u.m),
        )
        cartesian_cords = np.array(
            [coord.value for coord in self.cartesian_cords]
        )
        return N_fast(a, b, c, cartesian_cords)

    @property
    def tangential_vecs(self):
        """Returns orthonormal vectors tangential to the ellipsoid at the given location."""
        N = self.N
        return tangential_vecs_fast(N)

    @property
    def radius_of_curvature(self):
        """Radius of curvature of the meridian at the latitude of the given location."""
        _a, _c, _lat = (
            self._a.to_value(u.m),
            self._c.to_value(u.m),
            self._lat.to_value(u.rad),
        )
        return (
            radius_of_curvature_fast(_a, _c, _lat) * u.m
        )  # Need to convert units to u.rad and then take value because numpy expects angles in radians if unit is not given.

    def distance(self, px, py, pz):
        """
        Calculates the distance from an arbitrary point to the given location (Cartesian coordinates).

        Parameters
        ----------
        px : ~astropy.units.quantity.Quantity
            x-coordinate of the point
        py : ~astropy.units.quantity.Quantity
            y-coordinate of the point
        pz : ~astropy.units.quantity.Quantity
            z-coordinate of the point

        """
        px, py, pz = px.to_value(u.m), py.to_value(u.m), pz.to_value(u.m)
        cartesian_cords = np.array(
            [coord.value for coord in self.cartesian_cords]
        )
        return (
            distance_fast(cartesian_cords, px, py, pz) * u.m
        )  # body.R and body.R_polar has u.m as units

    def is_visible(self, px, py, pz):
        """
        Determines whether an object located at a given point is visible from the given location.
        Returns true if true, false otherwise.

        Parameters
        ----------
        px : ~astropy.units.quantity.Quantity
            x-coordinate of the point
        py : ~astropy.units.quantity.Quantity
            y-coordinate of the point
        pz : ~astropy.units.quantity.Quantity
            z-coordinate of the point

        """
        px, py, pz = px.to_value(u.m), py.to_value(u.m), pz.to_value(u.m)
        cartesian_cords = np.array(
            [coord.value for coord in self.cartesian_cords]
        )

        return is_visible_fast(cartesian_cords, px, py, pz, self.N)

    def cartesian_to_ellipsoidal(self, x, y, z):
        """
        Converts ellipsoidal coordinates to the Cartesian coordinate system for the given ellipsoid.

        Parameters
        ----------
        x : ~astropy.units.quantity.Quantity
            x coordinate
        y : ~astropy.units.quantity.Quantity
            y coordinate
        z : ~astropy.units.quantity.Quantity
            z coordinate

        """
        _a, _c = self._a.to_value(u.m), self._c.to_value(u.m)
        x, y, z = x.to_value(u.m), y.to_value(u.m), z.to_value(u.m)
        lat, lon, h = cartesian_to_ellipsoidal_fast(_a, _c, x, y, z)
        return lat * u.rad, lon * u.rad, h * u.m
