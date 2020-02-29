"""
Porkchops with poliastro
------------------------

Porkchops are also known as mission design curves since they show
different parameters used to design the ballistic trajectories for the
targetting problem such us:

-  Time of flight (TFL)
-  Launch energy (C3L)
-  Arrival velocity (VHP)

For the moment, poliastro is only capable of creating these mission
plots between ``poliastro.bodies`` objects. However, it is intended for
future versions to make it able for plotting porkchops between NEOs
also.

"""


######################################################################
# Basic modules
# ~~~~~~~~~~~~~
#
# For creating a porkchop plot with poliastro, we need to import the
# ``porkchop`` function from the ``poliastro.plotting.porkchop`` module.
# Also, two ``poliastro.bodies`` are necessary for computing the
# targetting problem associated. Finally by making use of ``time_range``,
# a very useful function available at ``poliastro.utils`` it is possible
# to define a span of launching and arrival dates for the problem.
#

import astropy.units as u
import matplotlib.pyplot as plt
from poliastro.bodies import Earth, Mars
from poliastro.plotting.porkchop import porkchop
from poliastro.util import time_range

launch_span = time_range("2005-04-30", end="2005-10-07")
arrival_span = time_range("2005-11-16", end="2006-12-21")


######################################################################
# Plot that porkchop!
# ~~~~~~~~~~~~~~~~~~~
#
# All that we must do is pass the two bodies, the two time spans and some
# extra plotting parameters realted to different information along the
# figure such us:
#
# -  If we want poliastro to plot time of flight lines: ``tfl=True/False``
# -  If we want poliastro to plot arrival velocity: ``vhp=True/False``
# -  The maximum value for C3 to be ploted:
#    ``max_c3=45 * u.km**2 / u.s**2`` (by default)
#

# We create the porkchop
dv_dpt, dv_arr, c3dpt, c3arr, tof = porkchop(Earth, Mars, launch_span, arrival_span)
plt.show()


######################################################################
# NASA's same porkchop
# ~~~~~~~~~~~~~~~~~~~~
#
# We can compare previous porkchop with the ones made by NASA for those
# years.
#
# .. figure:: images/porkchop_mars.png
#    :alt: Porkchop to Mars
#
#    Porkchop to Mars
#
