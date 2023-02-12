"""@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University.

This module complements constants defined in `astropy.constants` with
mean distances between sun-planet and planet-moon systems

Mean distances between the two bodies are required for CR3BP model CR3BP
model's goal is to approximate the behavior of a satellite in an ephemeris
model. Thus, it is reasonable to us low precision mean distances values as for
more comprehensive mission designs, the CR3BP trajectories are differentially
corrected in an ephemeris model with high precision mean distances.

Low precision mean distance between sun-planet were obtained from:

* NASA. (n.d.). Approximate positions of the planets. NASA., from
  https://ssd.jpl.nasa.gov/planets/approx_pos.html

Low precision mean distance between planet-moon were obtained from:

* NASA. (n.d.). Planetary satellite mean elements. NASA. from
  https://ssd.jpl.nasa.gov/sats/elem/sep.html

All computed at J2000, except "Pluto - Charon" -> Epoch: 2020-01-01.5 TDB

"""

from astropy import units as u
from astropy.constants import Constant

__all__ = [
    "mean_a_mercury",
    "mean_a_venus",
    "mean_a_earth",
    "mean_a_mars",
    "mean_a_jupiter",
    "mean_a_saturn",
    "mean_a_uranus",
    "mean_a_neptune",
    "mean_a_moon",
    "mean_a_phobos",
    "mean_a_deimos",
    "mean_a_europa",
    "mean_a_ganymede",
    "mean_a_enceladus",
    "mean_a_titan",
    "mean_a_titania",
    "mean_a_triton",
    "mean_a_charon",
]

mean_a_mercury = Constant(
    "mean_a_mercury",
    "Mercury's orbit mean semimajor axis around Sun",
    (0.38709927 * u.au.to(u.km)),
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

mean_a_venus = Constant(
    "mean_a_venus",
    "Venus's orbit mean semimajor axis around Sun",
    (0.72333566 * u.au.to(u.km)),
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

mean_a_earth = Constant(
    "mean_a_earth",
    "Earth's orbit mean semimajor axis around Sun",
    (1 * u.au.to(u.km)),
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

mean_a_mars = Constant(
    "mean_a_moon",
    "Mars's orbit mean semimajor axis around Sun",
    1.52371034 * u.au.to(u.km),
    "km",
    0,
    "NASA Approximate Position of Planets",
    system="si",
)

mean_a_jupiter = Constant(
    "mean_a_jupiter",
    "Jupiter's orbit mean semimajor axis around Sun",
    (5.20288700 * u.au.to(u.km)),
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

mean_a_saturn = Constant(
    "mean_a_saturn",
    "Saturn's orbit mean semimajor axis around Sun",
    (9.53667594 * u.au.to(u.km)),
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

mean_a_uranus = Constant(
    "mean_a_uranus",
    "Uranus' orbit mean semimajor axis around Sun",
    (19.18916464 * u.au.to(u.km)),
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

mean_a_neptune = Constant(
    "mean_a_neptune",
    "Neptune's orbit mean semimajor axis around Sun",
    (30.06992276 * u.au.to(u.km)),
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

mean_a_moon = Constant(
    "mean_a_moon",
    "Moon's orbit mean semimajor axis around Earth",
    384400,
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

mean_a_phobos = Constant(
    "mean_a_phobos",
    "Phobos's orbit mean semimajor axis around Mars",
    9400,
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

mean_a_deimos = Constant(
    "mean_a_deimos",
    "Deimos's orbit mean semimajor axis around Mars",
    23500,
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

mean_a_europa = Constant(
    "mean_a_europa",
    "Europa's orbit mean semimajor axis around Jupiter",
    671100,
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

mean_a_ganymede = Constant(
    "mean_a_ganymede",
    "Ganymede's orbit mean semimajor axis around Juputer",
    1070400,
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

mean_a_enceladus = Constant(
    "mean_a_enceladus",
    "Enceladus's orbit mean semimajor axis around Saturn",
    238400,
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

mean_a_titan = Constant(
    "mean_a_titan",
    "Titan's orbit mean semimajor axis around Saturn",
    1221900,
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

mean_a_titania = Constant(
    "mean_a_titania",
    "Titania's orbit mean semimajor axis around Uranus",
    436300,
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

mean_a_triton = Constant(
    "mean_a_triton",
    "Triton's orbit mean semimajor axis around Neptune",
    354800,
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)

mean_a_charon = Constant(
    "mean_a_charon",
    "Charon's orbit mean semimajor axis around Pluto",
    19600,
    "km",
    0,
    "NASA Planetary Satellite Mean Elements",
    system="si",
)
