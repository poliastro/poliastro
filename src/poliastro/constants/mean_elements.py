"""
@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University

This module complements constants defined in `astropy.constants` with
mean distances between sun-planet and planet-moon systems

Mean distances between the two bodies are required for CR3BP model
CR3BP model's goal is to approximate the behavior of a satellite in an
ephemeris model. Thus, it is reasonable to us low precision mean distances values as
for more comprehensive mission designs, the CR3BP trajectories are differentially
corrected in an ephemeris model with high precision mean distances.

Low precision mean distance between sun-planet were obtained from:
* NASA. (n.d.). Approximate positions of the planets. NASA.,
  from https://ssd.jpl.nasa.gov/planets/approx_pos.html

Low precision mean distance between planet-moon were obtained from:
* NASA. (n.d.). Planetary satellite mean elements. NASA.
  from https://ssd.jpl.nasa.gov/sats/elem/sep.html

All computed at J2000, except "Pluto - Charon" -> Epoch: 2020-01-01.5 TDB

"""

from astropy import units as u
from astropy.constants import Constant

__all__ = [
    "sunmercury",
    "sunvenus",
    "sunearth",
    "sunmars",
    "sunjupiter",
    "sunsaturn",
    "sunuranus",
    "sunneptune",
    "earthmoon",
    "marsphobos",
    "marsdeimos",
    "jupitereuropa",
    "jupiterganymede",
    "saturnenceladus",
    "saturntitan",
    "uranustitania",
    "neptunetriton",
    "plutocharon",
]

sunmercury = Constant(
    "sunmercury",
    "Mercury's orbit semimajor axis around Sun",
    (0.38709927 * u.au.to(u.km)),
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

sunvenus = Constant(
    "sunvenus",
    "Venus's orbit semimajor axis around Sun",
    (0.72333566 * u.au.to(u.km)),
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

sunearth = Constant(
    "sunearth",
    "Earth's orbit semimajor axis around Sun",
    (1 * u.au.to(u.km)),
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

sunmars = Constant(
    "earthmoon",
    "Mars's orbit semimajor axis around Sun",
    1.52371034 * u.au.to(u.km),
    "km",
    0,
    "NASA Approximate Position of Planets",
    system="si",
)

sunjupiter = Constant(
    "sunjupiter",
    "Jupiter's orbit semimajor axis around Sun",
    (5.20288700 * u.au.to(u.km)),
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

sunsaturn = Constant(
    "sunsaturn",
    "Saturn's orbit semimajor axis around Sun",
    (9.53667594 * u.au.to(u.km)),
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

sunuranus = Constant(
    "sunuranus",
    "Saturn's orbit semimajor axis around Sun",
    (19.18916464 * u.au.to(u.km)),
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

sunneptune = Constant(
    "sunneptune",
    "Neptune's orbit semimajor axis around Sun",
    (30.06992276 * u.au.to(u.km)),
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

earthmoon = Constant(
    "earthmoon",
    "Moon's orbit semimajor axis around Earth",
    384400,
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

marsphobos = Constant(
    "marsphobos",
    "Phobos's orbit semimajor axis around Mars",
    9400,
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

marsdeimos = Constant(
    "marsdeimos",
    "Deimos's orbit semimajor axis around Mars",
    23500,
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

jupitereuropa = Constant(
    "jupitereuropa",
    "Europa's orbit semimajor axis around Jupiter",
    671100,
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

jupiterganymede = Constant(
    "jupiterganymede",
    "Ganymede's orbit semimajor axis around Juputer",
    1070400,
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

saturnenceladus = Constant(
    "saturnenceladus",
    "Enceladus's orbit semimajor axis around Saturn",
    238400,
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

saturntitan = Constant(
    "saturntitan",
    "Titan's orbit semimajor axis around Saturn",
    1221900,
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

uranustitania = Constant(
    "uranustitania",
    "Titania's orbit semimajor axis around Uranus",
    436300,
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

neptunetriton = Constant(
    "neptunetriton",
    "Triton's orbit semimajor axis around Neptune",
    354800,
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

plutocharon = Constant(
    "plutocharon",
    "Charon's orbit semimajor axis around Pluto",
    19600,
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)
