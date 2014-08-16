# coding: utf-8
"""Example data.

"""

from astropy import time
from astropy import units as u
u.one = u.dimensionless_unscaled  # astropy #1980

from poliastro.bodies import Earth
from poliastro.twobody import State

# Taken from Plyades (c) 2012 Helge Eichhorn (MIT License)
iss = State.from_vectors(Earth,
                         [8.59072560e2, -4.13720368e3, 5.29556871e3] * u.km,
                         [7.37289205, 2.08223573, 4.39999794e-1] * u.km / u.s,
                         time.Time("2013-03-18 12:00", scale='utc'))

molniya = State.from_elements(Earth,
                              26600 * u.km, 0.75 * u.one, 63.4 * u.deg,
                              0 * u.deg, 270 * u.deg, 80 * u.deg)
