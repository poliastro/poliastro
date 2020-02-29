"""
Visualizing the SpaceX Tesla Roadster trip to Mars
==================================================

"""

from astropy.time import Time
from astropy import units as u
import matplotlib.pyplot as plt
import plotly

from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit
from poliastro.plotting import OrbitPlotter3D


EPOCH = Time("2018-02-18 12:00:00", scale="tdb")

roadster = Orbit.from_horizons("SpaceX Roadster", Sun, epoch=EPOCH, id_type="majorbody")
roadster

from poliastro.plotting.misc import plot_solar_system

frame = plot_solar_system(outer=False, epoch=EPOCH)
frame.plot(roadster, label="SpaceX Roadster", color="black");
plt.show()

######################################################################
# We will now make a three-dimensional plot
#

frame = OrbitPlotter3D()

frame.plot(Orbit.from_body_ephem(Earth, EPOCH), label=Earth)
frame.plot(Orbit.from_body_ephem(Mars, EPOCH), label=Mars)
frame.plot(roadster, label="SpaceX Roadster", color="black")

fig = frame.set_view(30 * u.deg, -100 * u.deg, 2 * u.km)
plotly.io.show(fig)

