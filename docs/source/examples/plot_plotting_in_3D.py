"""
Plotting in 3D
==============

"""

import numpy as np

from poliastro.examples import *
from poliastro.plotting import *

import plotly

######################################################################
# We start by plotting the Churi's orbit without the need of creating
# an instance of OrbitPlotter3D()
#

fig = churi.plot(interactive=True, use_3d=True)
plotly.io.show(fig)

######################################################################
# However, if we want to add more orbits to the frame, then we need
# an instance of the three-dimensional plotter utility
#

frame = OrbitPlotter3D()
frame.plot(churi)
fig = frame.plot(Orbit.from_body_ephem(Earth))
plotly.io.show(fig)


######################################################################
# Let us now plot for a Molniya orbit!
#

frame = OrbitPlotter3D()
frame.plot(molniya)
plotly.io.show(fig)


######################################################################
# What about adding the ISS to better compare both orbits?
#

frame = OrbitPlotter3D()

frame.plot(molniya)
fig = frame.plot(iss)
plotly.io.show(fig)


######################################################################
# We will now plot Eros' and Earth's orbit 
#

eros = Orbit.from_sbdb("eros")

frame = OrbitPlotter3D()
frame.plot(Orbit.from_body_ephem(Earth), label=Earth)
fig = frame.plot(eros, label="eros")
plotly.io.show(fig)


######################################################################
# We can also solve for a range of positions and velocities for the
# Earth at an specific period of time. When using the astropy
# get_body_barycentric_posvel() we need to set the attractor of the frame
# figure. 
#

from astropy.coordinates import get_body_barycentric_posvel
from poliastro.util import time_range

date_launch = time.Time("2011-11-26 15:02", scale="utc")
date_arrival = time.Time("2012-08-06 05:17", scale="utc")

rr_earth, _ = get_body_barycentric_posvel(
    "earth", time_range(date_launch, end=date_arrival, periods=50)
)

frame = OrbitPlotter3D()
frame.set_attractor(Sun)
frame.plot(Orbit.from_body_ephem(Earth), label=Earth)
fig = frame.plot_trajectory(rr_earth, label=Earth)
plotly.io.show(fig)


######################################################################
# However, if we first plot Eros, since an attractor has been set
# due to its Orbit class type, it is not necessary to specify the
# attractor of the plot anymore.

frame = OrbitPlotter3D()
frame.plot(eros, label="eros")
fig = frame.plot_trajectory(rr_earth, label=Earth)
plotly.io.show(fig)


