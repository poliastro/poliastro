import os
os.chdir("/home/arnau/Desktop/poliastro/src")
from astropy import units as u
from poliastro.plotting import *
from poliastro.twobody import Orbit
from poliastro.bodies import Earth
ss_i = Orbit.circular(Earth, alt=700 * u.km)
op = OrbitPlotter()
op.animate(ss_i)
plt.show()
