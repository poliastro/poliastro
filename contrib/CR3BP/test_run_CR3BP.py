"""
********************************************************************
      Test file for implementation check of CR3BP library.
********************************************************************

Last update: 21/01/2022

Description
-----------
Contains a few sample orbit propagations to test the CR3BP library.

The orbits currently found in test file include:
    - L2 southern NRHO (9:2 NRHO of Lunar Gateway Station)
    - Distant Retrograde Orbit (DRO)
    - Butterfly Orbit
    - L2 Vertical Orbit
"""

# Testing CR3BP implementation

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from CR3BP import getChar_CR3BP, propagate, propagateSTM

from poliastro.bodies import Earth, Moon

# Earth-Moon system properties
k1 = Earth.k.to(u.km**3 / u.s**2).value
k2 = Moon.k.to(u.km**3 / u.s**2).value
r12 = 384747.99198  # Earth-Moon distance

# Compute CR3BP characterisitic values
mu, kstr, lstr, tstr, vstr, nstr = getChar_CR3BP(k1, k2, r12)


# -- Lunar Gateway Station Orbit - 9:2 NRHO

"""
The orbit is a Near-Rectilinear Halo Orbit (NRHO) around the L2 Lagragian
point of the Earth-Moon system. The orbit presented here is a southern
sub-family of the L2-NRHO. This orbit is 9:2 resonant orbit currenly set
as the candidate orbit for the Lunar Gateway Station (LOP-G). Its called
9:2 resonant since a spacecraft would complete 9 orbits in the NRHO for
every 2 lunar month (slightly different from lunar orbit period).

The exact orbital elements presented here are from the auther's simulations.
The orbit states were obtained starting form guess solutions given in various
references. A few are provided below:

Ref: White Paper: Gateway Destination Orbit Model: A Continuous 15 Year NRHO
    Reference Trajectory - NASA, 2019
Ref: Strategies for Low-Thrust Transfer Design Based on Direct Collocation
    Techniques - Park, Howell and Folta

The NRHO are subfamily of the Halo orbits. The 'Near-Rectilinear' term comes
from the very elongated state of the orbit considering a regular Halo. Halo
orbits occur in all three co-linear equilibrum points L1,L2 and L3. They occur
in a pair of variants (nothern and southern) due to symmetry of CR3BP.
"""

# 9:2 L2 souther NRHO orbit
r0 = np.array([[1.021881345465263, 0, -0.182000000000000]])
v0 = np.array([0, -0.102950816739606, 0])
tf = 1.509263667286943

# number of points to plot
Nplt = 300
tofs = np.linspace(0, tf, Nplt)

# propagate the base trajectory
rf, vf = propagate(mu, r0, v0, tofs, rtol=1e-11)

# ploting orbit
rf = np.array(rf)

fig = plt.figure()
ax = plt.axes(projection="3d")
ax.set_box_aspect(
    (np.ptp(rf[:, 0]), np.ptp(rf[:, 1]), np.ptp(rf[:, 2]))
)  # aspect ratio is 1:1:1 in data space
# ploting the moon
ax.plot3D(1 - mu, 0, 0, "ok")
ax.set_title("L2 Southern NRHO")
ax.set_xlabel("x-axis [nd]")
ax.set_ylabel("y-axis [nd]")
ax.set_zlabel("z-axis [nd]")

ax.plot3D(rf[:, 0], rf[:, 1], rf[:, 2], "b")
plt.show()


"""
All other orbits in this section are computed from guess solutions available
in Grebow's Master and PhD thesis. He lists a quite detailed set of methods
to compute most of the major periodic orbits I have presented here. All of
them use differntial correction methods which are not yet implemented in this
library.

Ref: GENERATING PERIODIC ORBITS IN THE CIRCULAR RESTRICTED THREEBODY PROBLEM
    WITH APPLICATIONS TO LUNAR SOUTH POLE COVERAGE
    - D.Grebow 2006 (Master thesis)
Ref: TRAJECTORY DESIGN IN THE EARTH-MOON SYSTEM
    AND LUNAR SOUTH POLE COVERAGE
    - D.Grebow 2010 (PhD desertation)
"""


# -- DRO orbit

# DRO orbit states

r0 = np.array([0.783390492345344, 0, 0])
v0 = np.array([0, 0.548464515316651, 0])
tf = 3.63052604667440

# number of points to plot
Nplt = 300
tofs = np.linspace(0, tf, Nplt)

# propagate the base trajectory
rf, vf = propagate(mu, r0, v0, tofs, rtol=1e-11)


# ploting orbit
rf = np.array(rf)

fig = plt.figure()
ax = plt.axes(projection="3d")
ax.set_box_aspect(
    (np.ptp(rf[:, 0]), np.ptp(rf[:, 1]), np.ptp(rf[:, 2]))
)  # aspect ratio is 1:1:1 in data space
# ploting the moon
ax.plot3D(1 - mu, 0, 0, "ok")
ax.set_title("Distant Restrograde orbit (DRO)")
ax.set_xlabel("x-axis [nd]")
ax.set_ylabel("y-axis [nd]")
ax.set_zlabel("z-axis [nd]")

ax.plot3D(rf[:, 0], rf[:, 1], rf[:, 2], "m")
plt.show()


# -- Butterfly orbit

# Butterfly orbit states

r0 = np.array([1.03599510774957, 0, 0.173944812752286])
v0 = np.array([0, -0.0798042160573269, 0])
tf = 2.78676904546834

# number of points to plot
Nplt = 300
tofs = np.linspace(0, tf, Nplt)

# propagate the base trajectory
rf, vf = propagate(mu, r0, v0, tofs, rtol=1e-11)

# ploting orbit
rf = np.array(rf)

fig = plt.figure()
ax = plt.axes(projection="3d")
ax.set_box_aspect(
    (np.ptp(rf[:, 0]), np.ptp(rf[:, 1]), np.ptp(rf[:, 2]))
)  # aspect ratio is 1:1:1 in data space
# ploting the moon
ax.plot3D(1 - mu, 0, 0, "ok")
ax.set_title("Butterfly orbit")
ax.set_xlabel("x-axis [nd]")
ax.set_ylabel("y-axis [nd]")
ax.set_zlabel("z-axis [nd]")

ax.plot3D(rf[:, 0], rf[:, 1], rf[:, 2], "r")
plt.show()


# -- Vertical orbit

# Vertical orbit states

r0 = np.array([0.504689989562366, 0, 0.836429774762193])
v0 = np.array([0, 0.552722840538063, 0])
tf = 6.18448756121754

# number of points to plot
Nplt = 300
tofs = np.linspace(0, tf, Nplt)

# propagate the base trajectory
rf, vf = propagate(mu, r0, v0, tofs, rtol=1e-11)

# ploting orbit
rf = np.array(rf)

fig = plt.figure()
ax = plt.axes(projection="3d")
ax.set_box_aspect(
    (np.ptp(rf[:, 0]), np.ptp(rf[:, 1]), np.ptp(rf[:, 2]))
)  # aspect ratio is 1:1:1 in data space
# ploting the moon
ax.plot3D(1 - mu, 0, 0, "ok")
ax.set_title("L2 Vertical orbit")
ax.set_xlabel("x-axis [nd]")
ax.set_ylabel("y-axis [nd]")
ax.set_zlabel("z-axis [nd]")

ax.plot3D(rf[:, 0], rf[:, 1], rf[:, 2], "g")
plt.show()


# -- Propage STM

# propagate base trajectory with state-transition-matrix
STM0 = np.eye(6)
rf, vf, STM = propagateSTM(mu, r0, v0, STM0, tofs, rtol=1e-11)

# STM is a matrix of partial derivatives which are used in Newton-Raphson
# methods for trajectory design
