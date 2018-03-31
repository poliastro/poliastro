import os
os.chdir("/home/arnau/Desktop/poliastro/src")
from astropy import units as u
from poliastro.plotting import *
from poliastro.twobody import Orbit
from poliastro.bodies import Earth
ss_i = Orbit.circular(Earth, alt=700 * u.km)
ss_f = Orbit.circular(Earth, alt=1000 * u.km)
ss_n = Orbit.circular(Earth, alt=1500 * u.km)
op = OrbitPlotter()
op.plot(ss_i)
op.animate(ss_n)
op.animate(ss_f)
