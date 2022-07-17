"""
@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University
        dhruvj9922@gmail.com

Obj: To compute family of L4 Axial Orbit
    Single Shooter Variabl Time Setup
    1. Continue in 'vy' + Periodicity (sym_period_targ=1)
    2. PALC + Periodicity (sym_period_targ=1) {Phase Constraint is added}

Initial Condition obtained from:
sD. Grebow, "Generating Periodic Orbits in the Circular Restricted Three-Body Problem with Applications to Lunar South Pole Coverage," M.S., May 2006.
"""

import numpy as np
import plotly.graph_objs as go
from cr3bp_char_quant import sys_chars
from cr3bp_lib_calc import lib_pt_loc
from cr3bp_po_fam_continuation import periodic_orbit_fam_continuation
from cr3bp_po_plot_orbits import plot_orbits

sys_p1p2 = sys_chars('Earth','Moon')
mu = sys_p1p2.mu
lib_loc = lib_pt_loc(sys_p1p2)
li = lib_loc[3, :]  # 0 for L1 and  1 for L2


# From D. Grebow
ig = np.array([0.9522, -0.0146, 0.1, -0.1089, 0.0426, 0.4599])
tf_guess = 2.1684

orb_fam_obj = periodic_orbit_fam_continuation(sys_p1p2, ig,tf=tf_guess)

free_vars = ["x", "y", "z", "vx", "vy", "vz", "t"]
constraints = [
    "x",
    "y",
    "z",
    "vx",
    "vy",
]  # vz not constrained as JC implicitly constaints 5 of the 6 states

# Target Axial orbit using Single Shooter Variable Time setup
#        Exploits Periodicity(sym_perioid_targ set to 1)
#        Continue in 'vy' using Natural Paramter Continuaton
orb_fam_obj.npc_po_fam(free_vars, constraints,sym_period_targ=1, Nmax=10, 
                    step_size= 1e-3, num_fam_members=2, param_continue="vy", line_search=True)

'''
PALC
'''
orb_fam_obj.palc_po_fam(free_vars, constraints,sym_period_targ=1, Nmax=10, 
                    step_size= 7e-2, num_fam_members=100, line_search=True)

"""
Plot family
"""
# if targeted_po_char != None:
colourby = orb_fam_obj.targeted_po_char['jc']
colourmap='plasma'
cb_label = 'JC'
title = 'CR3BP: Earth-Moon L4 Northern Axial family: DJ'
data_trace = []
# Add L4
data_trace.append(go.Scatter3d(x=[li[0]], y=[0], z=[0], marker=dict(
            color='red',
            size=2)))
# Add Moon
data_trace.append(go.Scatter3d(x=[1-mu], y=[0], z=[0], marker=dict(
            color='grey',
            size=7)))
# Add Earth
data_trace.append(go.Scatter3d(x=[-mu], y=[0], z=[0], marker=dict(
            color='blue',
            size=10)))


plot_orbits(mu,orb_fam_obj.targeted_po_fam,colourby, cb_label, title=title,data_trace=data_trace, save=False)