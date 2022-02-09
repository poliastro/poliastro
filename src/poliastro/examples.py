"""Example data.

"""
from astropy import time, units as u

from poliastro.bodies import Earth, Sun
from poliastro.twobody import Orbit

iss = Orbit.from_vectors(
    Earth,
    [8.59072560e2, -4.13720368e3, 5.29556871e3] * u.km,
    [7.37289205, 2.08223573, 4.39999794e-1] * u.km / u.s,
    time.Time("2013-03-18 12:00", scale="utc"),
)
"""ISS orbit example

Taken from Plyades (c) 2012 Helge Eichhorn (MIT License)

"""

molniya = Orbit.from_classical(
    attractor=Earth,
    a=26600 * u.km,
    ecc=0.75 * u.one,
    inc=63.4 * u.deg,
    raan=0 * u.deg,
    argp=270 * u.deg,
    nu=80 * u.deg,
)
"""Molniya orbit example"""

_r_a = Earth.R + 35950 * u.km
_r_p = Earth.R + 250 * u.km
_a = (_r_a + _r_p) / 2
soyuz_gto = Orbit.from_classical(
    attractor=Earth,
    a=_a,
    ecc=_r_a / _a - 1,
    inc=6 * u.deg,
    raan=188.5 * u.deg,
    argp=178 * u.deg,
    nu=0 * u.deg,
)
"""Soyuz geostationary transfer orbit (GTO) example

Taken from Soyuz User's Manual, issue 2 revision 0

"""

churi = Orbit.from_classical(
    attractor=Sun,
    a=3.46250 * u.AU,
    ecc=0.64 * u.one,
    inc=7.04 * u.deg,
    raan=50.1350 * u.deg,
    argp=12.8007 * u.deg,
    nu=63.89 * u.deg,
    epoch=time.Time("2015-11-05 12:00", scale="utc"),
)
"""Comet 67P/Churyumov–Gerasimenko orbit example"""

debris_orbits_list = []

debris_orbits_list.append(Orbit.from_classical(
    attractor=Earth,
    a=6828.193338988509 * u.km,
    ecc=0.0062140354170737815 * u.one,
    inc=82.69387440482602 * u.deg,
    raan=37.33894561668519 * u.deg,
    argp=200.62393574484153 * u.deg,
    nu=-117.55203086408737 * u.deg
))

debris_orbits_list.append(Orbit.from_classical(
    attractor=Earth,
    a=6821.922877133498 * u.km,
    ecc=0.0023241234515223646 * u.one,
    inc=82.65684766470754 * u.deg,
    raan=36.3401421924121 * u.deg,
    argp=125.29597430617513 * u.deg,
    nu=-151.64963315597913 * u.deg
))

debris_orbits_list.append(Orbit.from_classical(
    attractor=Earth,
    a=6836.825166360441 * u.km,
    ecc=0.004635589373624103 * u.one,
    inc=82.69764910622918 * u.deg,
    raan=36.757861621556614 * u.deg,
    argp=44.219092511353594 * u.deg,
    nu=133.63349740950568 * u.deg
))

debris_orbits_list.append(Orbit.from_classical(
    attractor=Earth,
    a=7117.414488497581 * u.km,
    ecc=0.04540691741651712 * u.one,
    inc=83.07451144780156 * u.deg,
    raan=52.87995597314799 * u.deg,
    argp=190.63045916106168 * u.deg,
    nu=41.306044841636634 * u.deg
))

debris_orbits_list.append(Orbit.from_classical(
    attractor=Earth,
    a=6917.149445506758 * u.km,
    ecc=0.014811336719791061 * u.one,
    inc=82.6218902660939 * u.deg,
    raan=39.58175296436131 * u.deg,
    argp=106.71561062464224 * u.deg,
    nu=-62.66454424955413 * u.deg
))

debris_orbits_list.append(Orbit.from_classical(
    attractor=Earth,
    a=6910.923895002164 * u.km,
    ecc=0.016808366674815105 * u.one,
    inc=82.35369942440258 * u.deg,
    raan=35.60505049154483 * u.deg,
    argp=184.42211913686066 * u.deg,
    nu=45.95875318418421 * u.deg
))

debris_orbits_list.append(Orbit.from_classical(
    attractor=Earth,
    a=6677.880585440615 * u.km,
    ecc=0.0014802577911015675 * u.one,
    inc=82.17121703030627 * u.deg,
    raan=25.699484007134643 * u.deg,
    argp=204.23415215165576 * u.deg,
    nu=-24.40253410856961 * u.deg
))

debris_orbits_list.append(Orbit.from_classical(
    attractor=Earth,
    a=6961.126175813936 * u.km,
    ecc=0.025433234625720075 * u.one,
    inc=82.56625734212793 * u.deg,
    raan=40.073969157354824 * u.deg,
    argp=188.205744852877 * u.deg,
    nu=152.76011672297756 * u.deg
))

debris_orbits_list.append(Orbit.from_classical(
    attractor=Earth,
    a=6883.993995695628 * u.km,
    ecc=0.008022496129527943 * u.one,
    inc=83.2763699331857 * u.deg,
    raan=47.55215625476586 * u.deg,
    argp=114.15854367766484 * u.deg,
    nu=-31.793479924939778 * u.deg
))

debris_orbits_list.append(Orbit.from_classical(
    attractor=Earth,
    a=6798.16535966413 * u.km,
    ecc=0.006497072501655168 * u.one,
    inc=82.61152388022165 * u.deg,
    raan=34.66534775192231 * u.deg,
    argp=349.72219499585407 * u.deg,
    nu=90.16499314429353 * u.deg
))

debris_orbits_list.append(Orbit.from_classical(
    attractor=Earth,
    a=6854.939735394475 * u.km,
    ecc=0.006324699282668879 * u.one,
    inc=82.67294675705102 * u.deg,
    raan=37.411303678153935 * u.deg,
    argp=69.71516857133007 * u.deg,
    nu=-29.51423158409656 * u.deg
))
debris_orbits_list.append(Orbit.from_classical(
    attractor=Earth,
    a=6884.147529089194 * u.km,
    ecc=0.01552728578907267 * u.one,
    inc=82.4114164903979 * u.deg,
    raan=35.234082427651664 * u.deg,
    argp=186.60344739193755 * u.deg,
    nu=-146.9445890288543 * u.deg
))
debris_orbits_list.append(Orbit.from_classical(
    attractor=Earth,
    a=6846.946980540078 * u.km,
    ecc=0.0029371405696242913 * u.one,
    inc=82.6314212152875 * u.deg,
    raan=36.88448947562918 * u.deg,
    argp=134.53438085198738 * u.deg,
    nu=-69.91109773157386 * u.deg
))
debris_orbits_list.append(Orbit.from_classical(
    attractor=Earth,
    a=6914.901161035591 * u.km,
    ecc=0.010664583104155329 * u.one,
    inc=82.31703307692484 * u.deg,
    raan=34.654644407052835 * u.deg,
    argp=144.83292617925179 * u.deg,
    nu=-133.54025144695484 * u.deg
))
debris_orbits_list.append(Orbit.from_classical(
    attractor=Earth,
    a=7040.971575280624 * u.km,
    ecc=0.0333018067425175 * u.one,
    inc=82.50417227979605 * u.deg,
    raan=44.11015739081946 * u.deg,
    argp=133.5425169343891 * u.deg,
    nu=-42.74160359135228 * u.deg
))
"""Orbit List of COSMOS 1408 Debris Orbits Example

COSMOS 1408 Debris Data Taken from https://celestrak.com/NORAD/

"""
