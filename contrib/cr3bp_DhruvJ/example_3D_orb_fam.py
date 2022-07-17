"""
Created on Sat Feb 26 23:09:08 2022

@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University
        dhruvj9922@gmail.com

Obj: Test Periodic Orbit Single Shooter

Initial Condition obtained from:
D. Grebow, "Generating Periodic Orbits in the Circular Restricted Three-Body Problem with Applications to Lunar South Pole Coverage," M.S., May 2006.
"""

import copy

import matplotlib.pyplot as plt
import numpy as np
from cr3bp_char_quant import sys_chars
from cr3bp_lib_calc import lib_pt_loc
from cr3bp_PO_master import periodic_orbit

sys_p1p2 = sys_chars("Earth", "Moon")
mu = sys_p1p2.mu
lib_loc = lib_pt_loc(sys_p1p2)
li = lib_loc[1, :]  # 0 for L1 and  1 for L2

ig = np.array(
    [1.021881345465263, 0, 0.182000000000000, 0, -0.102950816739606, 0]
)
tf_guess = 1.509263667286943

free_vars = ["x", "vy", "t"]
constraints = ["y", "vx", "vz"]

orbit_results = []

po = periodic_orbit(sys_p1p2, ig)
po.tf = tf_guess

for i in range(20):
    results, iterflag = po.single_shooter(free_vars, constraints)

    orbit_results.append(results)

    po.ic = copy.copy(results["states"][0, :])
    print(po.ic)
    po.ic[2] += 0.001


# results, iterflag = po_single_shooter_cr3bp(mu, ig, tf_guess, free_vars, constraints)
# orbit_results.append(results)

plt.figure(1)
ax = plt.axes(projection="3d")
ax.set_title("EM, L1 Orbit Family, tol = 1e-12")
for i in range(len(orbit_results)):
    ax.plot3D(
        orbit_results[i]["states"][:, 0],
        orbit_results[i]["states"][:, 1],
        orbit_results[i]["states"][:, 2],
    )
plt.plot(li[0], li[1], "ro", label="L1")
ax.scatter(li[0], li[1], li[2], color="red")
ax.set_box_aspect(
    [ub - lb for lb, ub in (getattr(ax, f"get_{a}lim")() for a in "xyz")]
)
ax.set_ylabel("y [nd]")
ax.set_xlabel("x [nd]")
ax.set_zlabel("z [nd]")
plt.show()
