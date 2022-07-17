"""
@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University
        dhruvj9922@gmail.com

Obj: To compute family of L3 Axial Orbit
    Single Shooter Variabl Time Setup
    1. Continue in 'vz' + X-axis symmetry => targets Period/2 states
    2. PALC+ X-axis symmetry => targets Period/2 states

Initial Condition obtained from:
D. Grebow, "Generating Periodic Orbits in the Circular Restricted Three-Body Problem with Applications to Lunar South Pole Coverage," M.S., May 2006.
"""

import numpy as np
import plotly.graph_objs as go
from cr3bp_char_quant import sys_chars
from cr3bp_lib_calc import lib_pt_loc
from cr3bp_po_fam_continuation import periodic_orbit_fam_continuation
from cr3bp_po_plot_orbits import plot_orbits

sys_p1p2 = sys_chars("Earth", "Moon")
mu = sys_p1p2.mu
lib_loc = lib_pt_loc(sys_p1p2)
li = lib_loc[2, :]  # 0 for L1 and  1 for L2


# From D. Grebow
ig = np.array([-1.8963, 0, 0, 0, 1.6715, 1e-16])
tf_guess = 6.2616

orb_fam_obj = periodic_orbit_fam_continuation(sys_p1p2, ig, tf=tf_guess)

free_vars = ["x", "vz", "vy", "t"]
constraints = ["y", "z", "vx"]

# Target Axial orbit using Single Shooter Variable Time setup
#        Exploits X-axis symmetry(sym_perioid_targ set to 1/2)
#        Continue in 'vz' using Natural Paramter Continuaton
orb_fam_obj.npc_po_fam(
    free_vars,
    constraints,
    sym_period_targ=1 / 2,
    Nmax=10,
    step_size=1e-5,
    num_fam_members=3,
    param_continue="vz",
    line_search=True,
)

"""
PALC
# """
# Target Axial orbit using Single Shooter Variable Time setup
#        Exploits X-axis symmetry(sym_perioid_targ set to 1/2)
orb_fam_obj.palc_po_fam(
    free_vars,
    constraints,
    sym_period_targ=1 / 2,
    Nmax=10,
    step_size=1e-2 * 5,
    num_fam_members=10,
    line_search=True,
)


"""
Plot family
"""
# if targeted_po_char != None:
colourby = orb_fam_obj.targeted_po_char["jc"]
colourmap = "plasma"
cb_label = "JC"
title = "EM_L3_Axial_family_PALC"
data_trace = []
# Add L2
data_trace.append(
    go.Scatter3d(x=[li[0]], y=[0], z=[0], marker=dict(color="red", size=2))
)
# Add Earth
data_trace.append(
    go.Scatter3d(x=[1 - mu], y=[0], z=[0], marker=dict(color="grey", size=7))
)

plot_orbits(
    mu,
    orb_fam_obj.targeted_po_fam,
    colourby,
    cb_label,
    title=title,
    data_trace=data_trace,
    save=False,
)
