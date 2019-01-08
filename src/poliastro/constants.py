"""Astronomical and physics constants.

This module complements constants defined in `astropy.constants`,
with gravitational paremeters and radii.

Note that `GM_jupiter` and `GM_neptune` are both referred to the whole planetary system gravitational parameter.

Unless otherwise specified, gravitational and mass parameters were obtained from:

* Luzum, Brian et al. “The IAU 2009 System of Astronomical Constants: The Report of the IAU Working Group on Numerical
  Standards for Fundamental Astronomy.” Celestial Mechanics and Dynamical Astronomy 110.4 (2011): 293–304.
  Crossref. Web. `DOI: 10.1007/s10569-011-9352-4`_

radii were obtained from:

* Archinal, B. A. et al. “Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009.”
  Celestial Mechanics and Dynamical Astronomy 109.2 (2010): 101–135. Crossref. Web. `DOI: 10.1007/s10569-010-9320-4`_

.. _`DOI: 10.1007/s10569-011-9352-4`: http://dx.doi.org/10.1007/s10569-011-9352-4
.. _`DOI: 10.1007/s10569-010-9320-4`: http://dx.doi.org/10.1007/s10569-010-9320-4

J2 for the Sun was obtained from:

* https://hal.archives-ouvertes.fr/hal-00433235/document (New values of gravitational moments J2 and J4 deduced
  from helioseismology, Redouane Mecheri et al)

"""
from astropy import time
from astropy.constants import Constant

# See for example USNO Circular 179
J2000 = time.Time("J2000", scale="tt")

GM_sun = Constant(
    "GM_sun",
    "Heliocentric gravitational constant",
    1.32712442099e20,
    "m3 / (s2)",
    0.0000000001e20,
    "IAU 2009 system of astronomical constants",
    system="si",
)

GM_earth = Constant(
    "GM_earth",
    "Geocentric gravitational constant",
    3.986004418e14,
    "m3 / (s2)",
    0.000000008e14,
    "IAU 2009 system of astronomical constants",
    system="si",
)

# Anderson, John D. et al. “The Mass, Gravity Field, and Ephemeris of Mercury.” Icarus 71.3 (1987): 337–349.
# Crossref. Web. DOI: 10.1016/0019-1035(87)90033-9
GM_mercury = Constant(
    "GM_mercury",
    "Mercury gravitational constant",
    2.203209e13,
    "m3 / (s2)",
    0.91,
    "IAU 2009 system of astronomical constants",
    system="si",
)

# Konopliv, A.S., W.B. Banerdt, and W.L. Sjogren. “Venus Gravity: 180th Degree and Order Model.”
# Icarus 139.1 (1999): 3–18. Crossref. Web. DOI: 10.1006/icar.1999.6086
GM_venus = Constant(
    "GM_venus",
    "Venus gravitational constant",
    3.24858592e14,
    "m3 / (s2)",
    0.006,
    "IAU 2009 system of astronomical constants",
    system="si",
)

# Konopliv, Alex S. et al. “A Global Solution for the Mars Static and Seasonal Gravity, Mars Orientation, Phobos and
# Deimos Masses, and Mars Ephemeris.” Icarus 182.1 (2006): 23–50.
# Crossref. Web. DOI: 10.1016/j.icarus.2005.12.025
GM_mars = Constant(
    "GM_mars",
    "Mars gravitational constant",
    4.282837440e13,
    "m3 / (s2)",
    0.00028,
    "IAU 2009 system of astronomical constants",
    system="si",
)

# Jacobson, R. A. et al. “A comprehensive orbit reconstruction for the galileo prime mission in the JS200 system.”
# The Journal of the Astronautical Sciences 48.4 (2000): 495–516.
# Crossref. Web.
GM_jupiter = Constant(
    "GM_jupiter",
    "Jovian system gravitational constant",
    1.2671276253e17,
    "m3 / (s2)",
    2.00,
    "IAU 2009 system of astronomical constants",
    system="si",
)

# Jacobson, R. A. et al. “The Gravity Field of the Saturnian System from Satellite Observations and Spacecraft
# Tracking Data.” The Astronomical Journal 132.6 (2006): 2520–2526.
# Crossref. Web. DOI: 10.1086/508812
GM_saturn = Constant(
    "GM_saturn",
    "Saturn gravitational constant",
    3.79312077e16,
    "m3 / (s2)",
    1.1,
    "IAU 2009 system of astronomical constants",
    system="si",
)

# Jacobson, R. A. et al. “The Masses of Uranus and Its Major Satellites from Voyager Tracking Data and Earth-Based
# Uranian Satellite Data.” The Astronomical Journal 103 (1992): 2068.
# Crossref. Web. DOI: 10.1086/116211
GM_uranus = Constant(
    "GM_uranus",
    "Uranus gravitational constant",
    5.7939393e15,
    "m3 / (s2)",
    13.0,
    "IAU 2009 system of astronomical constants",
    system="si",
)

# Jacobson, R. A. “THE ORBITS OF THE NEPTUNIAN SATELLITES AND THE ORIENTATION OF THE POLE OF NEPTUNE.”
# The Astronomical Journal 137.5 (2009): 4322–4329. Crossref. Web. DOI:
# 10.1088/0004-6256/137/5/4322
GM_neptune = Constant(
    "GM_neptune",
    "Neptunian system gravitational constant",
    6.836527100580397e15,
    "m3 / (s2)",
    10.0,
    "IAU 2009 system of astronomical constants",
    system="si",
)

# Tholen, David J. et al. “MASSES OF NIX AND HYDRA.” The Astronomical Journal 135.3 (2008): 777–784. Crossref. Web.
# DOI: 10.1088/0004-6256/135/3/777
GM_pluto = Constant(
    "GM_pluto",
    "Pluto gravitational constant",
    8.703e11,
    "m3 / (s2)",
    3.7,
    "IAU 2009 system of astronomical constants",
    system="si",
)

# Lemoine, Frank G. et al. “High-Degree Gravity Models from GRAIL Primary Mission Data.”
# Journal of Geophysical Research: Planets 118.8 (2013): 1676–1698.
# Crossref. Web. DOI: 10.1002/jgre.20118
GM_moon = Constant(
    "GM_moon",
    "Moon gravitational constant",
    4.90279981e12,
    "m3 / (s2)",
    0.00000774,
    "Journal of Geophysical Research: Planets 118.8 (2013)",
    system="si",
)


R_sun = Constant(
    "R_sun",
    "Sun equatorial radius",
    6.96000e8,
    "m",
    0,
    "IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009",
    system="si",
)

R_earth = Constant(
    "R_earth",
    "Earth equatorial radius",
    6.3781366e6,
    "m",
    0.0001,
    "IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009",
    system="si",
)

R_mercury = Constant(
    "R_mercury",
    "Mercury equatorial radius",
    2.4397e6,
    "m",
    1.0,
    "IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009",
    system="si",
)

R_venus = Constant(
    "R_venus",
    "Venus equatorial radius",
    6.0518e6,
    "m",
    1.0,
    "IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009",
    system="si",
)

R_mars = Constant(
    "R_mars",
    "Mars equatorial radius",
    3.39619e6,
    "m",
    0.1,
    "IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009",
    system="si",
)

R_jupiter = Constant(
    "R_jupiter",
    "Jupiter equatorial radius",
    7.1492e7,
    "m",
    4,
    "IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009",
    system="si",
)

R_saturn = Constant(
    "R_saturn",
    "Saturn equatorial radius",
    6.0268e7,
    "m",
    4,
    "IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009",
    system="si",
)

R_uranus = Constant(
    "R_uranus",
    "Uranus equatorial radius",
    2.5559e7,
    "m",
    4,
    "IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009",
    system="si",
)

R_neptune = Constant(
    "R_neptune",
    "Neptune equatorial radius",
    2.4764e7,
    "m",
    15,
    "IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009",
    system="si",
)

R_pluto = Constant(
    "R_pluto",
    "Pluto effective radius",
    1.195e6,
    "m",
    5,
    "IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009",
    system="si",
)

R_moon = Constant(
    "R_moon",
    "Moon equatorial radius",
    1.7374e6,
    "m",
    1,
    "IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009",
    system="si",
)

J2_sun = Constant(
    "J2_sun",
    "Sun J2 non-oblateness coefficient",
    2.20e-7,
    "",
    0.01e-7,
    "HAL archives",
    system="si",
)

J2_earth = Constant(
    "J2_earth",
    "Earth J2 non-oblateness coefficient",
    0.00108263,
    "",
    1,
    "HAL archives",
    system="si",
)

J3_earth = Constant(
    "J3_earth",
    "Earth J3 non-oblateness coefficient",
    -2.5326613168e-6,
    "",
    1,
    "HAL archives",
    system="si",
)

H0_earth = Constant(
    "H0_earth",
    "Earth H0 atmospheric scale height",
    8500,
    "m",
    1,
    "de Pater and Lissauer 2010",
    system="si",
)

rho0_earth = Constant(
    "rho0_earth",
    "Earth rho0 atmospheric density prefactor",
    1.3,
    "kg / (m3)",
    1,
    "de Pater and Lissauer 2010",
    system="si",
)

Wdivc_sun = Constant(
    "Wdivc_sun",
    "total radiation power of Sun divided by the speed of light",
    1.0203759306204136e14,
    "kg km / (s2)",
    1,
    "Howard Curtis",
    system="si",
)
