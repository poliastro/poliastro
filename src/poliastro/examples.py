# coding: utf-8
"""Example data.

"""
from astropy import time
from astropy import units as u

from poliastro.bodies import Sun, Earth
from poliastro.twobody import Orbit


# Taken from Plyades (c) 2012 Helge Eichhorn (MIT License)
iss = Orbit.from_vectors(
    Earth,
    [8.59072560e2, -4.13720368e3, 5.29556871e3] * u.km,
    [7.37289205, 2.08223573, 4.39999794e-1] * u.km / u.s,
    time.Time("2013-03-18 12:00", scale='utc')
)

molniya = Orbit.from_classical(
    Earth,
    26600 * u.km, 0.75 * u.one, 63.4 * u.deg,
    0 * u.deg, 270 * u.deg, 80 * u.deg
)

churi = Orbit.from_classical(
    Sun,
    3.46250 * u.AU, 0.64 * u.one, 7.04 * u.deg,
    50.1350 * u.deg, 12.8007 * u.deg, 63.89 * u.deg,
    time.Time("2015-11-05 12:00", scale='utc')
)
