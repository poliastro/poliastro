# Testing CR3BP implementation

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

from poliastro.bodies import Earth, Moon
from CR3BP import getChar_CR3BP, propagate, propagateSTM
from astropy import units as u

# Earth-Moon system properties
k1 = Earth.k.to(u.km**3/u.s**2).value
k2 = Moon.k.to(u.km**3/u.s**2).value
r12 = 384747.99198  # Earth-Moon distance

# Compute CR3BP characterisitic values
mu, kstr, lstr, tstr, vstr, nstr = getChar_CR3BP(k1,k2,r12)



# -- Lunar Gateway Station Orbit - 9:2 NRHO

# 9:2 L2 souther NRHO orbit
r0 = np.array([[1.021881345465263,0,-0.182000000000000]])
v0 = np.array([0,-0.102950816739606,0])
tf = 1.509263667286943

# number of points to plot
Nplt = 300
tofs = np.linspace(0,tf,Nplt)

# propagate the base trajectory
rf,vf = propagate(mu, r0, v0, tofs, rtol=1e-11)

# ploting orbit
rf = np.array(rf)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_box_aspect((np.ptp(rf[:,0]), np.ptp(rf[:,1]), np.ptp(rf[:,2])))  # aspect ratio is 1:1:1 in data space
# ploting the moon
ax.plot3D(1-mu,0,0,'ok')
ax.set_title("L2 Southern NRHO")
ax.set_xlabel("x-axis [nd]")
ax.set_ylabel("y-axis [nd]")
ax.set_zlabel("z-axis [nd]")

ax.plot3D(rf[:,0],rf[:,1],rf[:,2],'b')



# -- DRO orbit

# DRO orbit states
	
r0 = np.array([0.783390492345344, 0, 0])
v0 = np.array([0, 0.548464515316651, 0])
tf = 3.63052604667440

# number of points to plot
Nplt = 300
tofs = np.linspace(0,tf,Nplt)

# propagate the base trajectory
rf,vf = propagate(mu, r0, v0, tofs, rtol=1e-11)


# ploting orbit
rf = np.array(rf)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_box_aspect((np.ptp(rf[:,0]), np.ptp(rf[:,1]), np.ptp(rf[:,2])))  # aspect ratio is 1:1:1 in data space
# ploting the moon
ax.plot3D(1-mu,0,0,'ok')
ax.set_title("Distant Restrograde orbit (DRO)")
ax.set_xlabel("x-axis [nd]")
ax.set_ylabel("y-axis [nd]")
ax.set_zlabel("z-axis [nd]")

ax.plot3D(rf[:,0],rf[:,1],rf[:,2],'m')





# -- Butterfly orbit

# Butterfly orbit states
	
r0 = np.array([1.03599510774957, 0, 0.173944812752286])
v0 = np.array([0, -0.0798042160573269, 0])
tf = 2.78676904546834

# number of points to plot
Nplt = 300
tofs = np.linspace(0,tf,Nplt)

# propagate the base trajectory
rf,vf = propagate(mu, r0, v0, tofs, rtol=1e-11)

# ploting orbit
rf = np.array(rf)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_box_aspect((np.ptp(rf[:,0]), np.ptp(rf[:,1]), np.ptp(rf[:,2])))  # aspect ratio is 1:1:1 in data space
# ploting the moon
ax.plot3D(1-mu,0,0,'ok')
ax.set_title("Butterfly orbit")
ax.set_xlabel("x-axis [nd]")
ax.set_ylabel("y-axis [nd]")
ax.set_zlabel("z-axis [nd]")

ax.plot3D(rf[:,0],rf[:,1],rf[:,2],'m')



# -- Vertical orbit

# Vertical orbit states
	
r0 = np.array([0.504689989562366, 0, 0.836429774762193])
v0 = np.array([0, 0.552722840538063, 0])
tf = 6.18448756121754

# number of points to plot
Nplt = 300
tofs = np.linspace(0,tf,Nplt)

# propagate the base trajectory
rf,vf = propagate(mu, r0, v0, tofs, rtol=1e-11)

# ploting orbit
rf = np.array(rf)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_box_aspect((np.ptp(rf[:,0]), np.ptp(rf[:,1]), np.ptp(rf[:,2])))  # aspect ratio is 1:1:1 in data space
# ploting the moon
ax.plot3D(1-mu,0,0,'ok')
ax.set_title("L2 Vertical orbit")
ax.set_xlabel("x-axis [nd]")
ax.set_ylabel("y-axis [nd]")
ax.set_zlabel("z-axis [nd]")

ax.plot3D(rf[:,0],rf[:,1],rf[:,2],'m')



# -- Propage STM

# propagate base trajectory with state-transition-matrix
STM0 = np.eye(6)
rf,vf,STM = propagateSTM(mu, r0, v0, STM0, tofs, rtol=1e-11)

# STM is a matrix of partial derivatives which are used in Newton-Raphson
# methods for trajectory design
