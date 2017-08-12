"""Astronomical and physics constants in SI units.

This module complements constants defined in `astropy.constants`

Unless otherwise specified, gravitational and mass parameters were obtained from:

* Archinal, B. A. et al. “Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009.”
  Celestial Mechanics and Dynamical Astronomy 109.2 (2010): 101–135. Crossref. Web. `DOI: 10.1007/s10569-010-9320-4`_

and radii were obtained from:

* Luzum, Brian et al. “The IAU 2009 System of Astronomical Constants: The Report of the IAU Working Group on Numerical
  Standards for Fundamental Astronomy.” Celestial Mechanics and Dynamical Astronomy 110.4 (2011): 293–304.
  Crossref. Web. `DOI: 10.1007/s10569-011-9352-4`_


.. _`DOI: 10.1007/s10569-010-9320-4`: http://dx.doi.org/10.1007/s10569-010-9320-4
.. _`DOI: 10.1007/s10569-011-9352-4`: http://dx.doi.org/10.1007/s10569-011-9352-4

"""

from astropy.constants import Constant


G = Constant('G', 'Constant of gravitation', 6.67428e-11, 'm3 / (kg * s2)', 6.7e-15,
             'IAU 2009 system of astronomical constants', system='si')

GM_sun = Constant('GM_sun', 'Heliocentric gravitational constant', 1.32712442099e20, 'm3 / (s2)', 1.0e10,
                  'IAU 2009 system of astronomical constants', system='si')

GM_earth = Constant('GM_earth', 'Geocentric gravitational constant', 3.986004418e14, 'm3 / (s2)', 8.0e5,
                    'IAU 2009 system of astronomical constants', system='si')

# Anderson, John D. et al. “The Mass, Gravity Field, and Ephemeris of Mercury.” Icarus 71.3 (1987): 337–349.
# Crossref. Web. DOI: 10.1016/0019-1035(87)90033-9
GM_mercury = Constant('GM_mercury', 'Mercury gravitational constant', 22032.09, 'km3 / (s2)', 0.91,
                      'IAU 2009 system of astronomical constants')

# Konopliv, A.S., W.B. Banerdt, and W.L. Sjogren. “Venus Gravity: 180th Degree and Order Model.”
# Icarus 139.1 (1999): 3–18. Crossref. Web. DOI: 10.1006/icar.1999.6086
GM_venus = Constant('GM_venus', 'Venus gravitational constant', 324858.592, 'km3 / (s2)', 0.006,
                    'IAU 2009 system of astronomical constants')

# Konopliv, Alex S. et al. “A Global Solution for the Mars Static and Seasonal Gravity, Mars Orientation, Phobos and
# Deimos Masses, and Mars Ephemeris.” Icarus 182.1 (2006): 23–50. Crossref. Web. DOI: 10.1016/j.icarus.2005.12.025
GM_mars = Constant('GM_mars', 'Mars gravitational constant', 42828.37440, 'km3 / (s2)', 0.00028,
                   'IAU 2009 system of astronomical constants')

# Jacobson, R. A. et al. “A comprehensive orbit reconstruction for the galileo prime mission in the JS200 system.”
# The Journal of the Astronautical Sciences 48.4 (2000): 495–516. Crossref. Web.
GM_jupiter = Constant('GM_jupiter', 'Jovian system gravitational constant', 126712762.53, 'km3 / (s2)', 2.00,
                      'IAU 2009 system of astronomical constants')

# Jacobson, R. A. et al. “The Gravity Field of the Saturnian System from Satellite Observations and Spacecraft
# Tracking Data.” The Astronomical Journal 132.6 (2006): 2520–2526. Crossref. Web. DOI: 10.1086/508812
GM_saturn = Constant('GM_saturn', 'Saturn gravitational constant', 37931207.7, 'km3 / (s2)', 1.1,
                     'IAU 2009 system of astronomical constants')

# Jacobson, R. A. et al. “The Masses of Uranus and Its Major Satellites from Voyager Tracking Data and Earth-Based
# Uranian Satellite Data.” The Astronomical Journal 103 (1992): 2068. Crossref. Web. DOI: 10.1086/116211
GM_uranus = Constant('GM_uranus', 'Uranus gravitational constant', 5793939.3, 'km3 / (s2)', 13.0,
                     'IAU 2009 system of astronomical constants')

# Jacobson, R. A. “THE ORBITS OF THE NEPTUNIAN SATELLITES AND THE ORIENTATION OF THE POLE OF NEPTUNE.”
# The Astronomical Journal 137.5 (2009): 4322–4329. Crossref. Web. DOI: 10.1088/0004-6256/137/5/4322
GM_neptune = Constant('GM_neptune', 'Neptunian system gravitational constant', 6836527.100580397, 'km3 / (s2)', 10.0,
                      'IAU 2009 system of astronomical constants')

# Tholen, David J. et al. “MASSES OF NIX AND HYDRA.” The Astronomical Journal 135.3 (2008): 777–784. Crossref. Web.
# DOI: 10.1088/0004-6256/135/3/777
GM_pluto = Constant('GM_pluto', 'Pluto gravitational constant', 870.3, 'km3 / (s2)', 3.7,
                    'IAU 2009 system of astronomical constants')

# Lemoine, Frank G. et al. “High-Degree Gravity Models from GRAIL Primary Mission Data.”
# Journal of Geophysical Research: Planets 118.8 (2013): 1676–1698. Crossref. Web. DOI: 10.1002/jgre.20118
GM_moon = Constant('GM_moon', 'Moon gravitational constant', 4902.79981, 'km3 / (s2)', 7.74e-06,
                   'Journal of Geophysical Research: Planets 118.8 (2013)')


R_sun = Constant('R_sun', 'Sun equatorial radius', 696000, 'km', 0,
                 'IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009')

R_earth = Constant('R_earth', 'Earth equatorial radius', 6378.1366, 'km', 0.0001,
                   'IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009')

R_mercury = Constant('R_mercury', 'Mercury equatorial radius', 2439.7, 'km', 1.0,
                     'IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009')

R_venus = Constant('R_venus', 'Venus equatorial radius', 6051.8, 'km', 1.0,
                   'IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009')

R_mars = Constant('R_mars', 'Mars equatorial radius', 3396.19, 'km', 0.1,
                  'IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009')

R_jupiter = Constant('R_jupiter', 'Jupiter equatorial radius', 71492, 'km', 4,
                     'IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009')

R_saturn = Constant('R_saturn', 'Saturn equatorial radius', 60268, 'km', 4,
                    'IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009')

R_uranus = Constant('R_uranus', 'Uranus equatorial radius', 25559, 'km', 4,
                    'IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009')

R_neptune = Constant('R_neptune', 'Neptune equatorial radius', 24764, 'km', 15,
                     'IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009')

R_pluto = Constant('R_pluto', 'Pluto effective radius', 1195, 'km', 5,
                   'IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009')

R_moon = Constant('R_moon', 'Moon equatorial radius', 1737.4, 'km', 1,
                  'IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009')
