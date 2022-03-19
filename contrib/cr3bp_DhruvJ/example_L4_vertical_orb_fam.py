"""
@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University
        dhruvj9922@gmail.com

Obj: To compute family of L2 Vertical Orbit
    Single Shooter Variabl Time Setup
    1. Continue in 'x' + XZ plane symmetry and X-axis symmetry use => targets Period/4 states
    2. Continue in 'x' + Periodicity targeter => targets Period states
    1. Continue in 'jc' + XZ plane symmetry and X-axis symmetry use => targets Period/4 states

Initial Condition obtained from:
D. Grebow, "Generating Periodic Orbits in the Circular Restricted Three-Body Problem with Applications to Lunar South Pole Coverage," M.S., May 2006.
"""
import numpy as np
import plotly.graph_objs as go
from cr3bp_char_quant import bodies_char
from cr3bp_fam_continuation import npc_po_fam_cr3bp
from cr3bp_lib_JC_calc import lib_pt_loc
from cr3bp_plot_orbits import plot_orbits
from cr3bp_PO_targeter import po_single_shooter_cr3bp

mu, dist_e_m, tstar = bodies_char("Earth", "Moon")
lib_loc = lib_pt_loc(mu)
li = lib_loc[1, :]  # 0 for L1 and  1 for L2

# From D. Grebow
ig = np.array([1.0842, 0, 0, 0, -0.5417, 0.8415])
tf_guess = 6.1305

orbit_results = []

free_vars = ["x", "vy", "vz", "t"]
constraints = ["y", "vx", "vz"]

# Target Vertical orbit using Single Shooter Variable Time setup
#        Exploits XZ plane symmetry and X-axis symmetry(sym_perioid_targ set to 1/4)
#        Continue in 'x' using Natural Paramter Continuaton to compute 20 family members
targeted_po_fam, targeted_po_char = npc_po_fam_cr3bp(
    mu,
    po_single_shooter_cr3bp,
    ig,
    tf_guess,
    free_vars,
    constraints,
    sym_period_targ=1 / 4,
    JCd=None,
    step_size=-1e-4,
    num_fam_members=20,
    param_continue="x",
    line_search=True,
)

constraints = ["y", "x", "vz"]
# Target Vertical orbit using Single Shooter Variable Time setup
#        Exploits Periodcity(sym_perioid_targ set to 1)
#        Continue in 'x' using Natural Paramter Continuaton to compute 20 family members
initial_guess = targeted_po_fam[-1]["states"][0, :]
tf_guess = targeted_po_fam[-1]["t"][-1]
targeted_po_fam_updated, targeted_po_char_updated = npc_po_fam_cr3bp(
    mu,
    po_single_shooter_cr3bp,
    ig,
    tf_guess,
    free_vars,
    constraints,
    sym_period_targ=1,
    JCd=None,
    step_size=-1e-4,
    num_fam_members=20,
    param_continue="x",
    line_search=True,
)
targeted_po_fam.extend(targeted_po_fam_updated)
for keys in targeted_po_char_updated.keys():
    targeted_po_char[keys].extend(targeted_po_char_updated[keys])

constraints = ["y", "vx", "vz"]
# Target Vertical orbit using Single Shooter Variable Time setup
#        Exploits Periodcity(sym_perioid_targ set to 1)
#        Continue in 'x' using Natural Paramter Continuaton to compute 20 family members
initial_guess = targeted_po_fam[-1]["states"][0, :]
tf_guess = targeted_po_fam[-1]["t"][-1]
targeted_po_fam_updated, targeted_po_char_updated = npc_po_fam_cr3bp(
    mu,
    po_single_shooter_cr3bp,
    ig,
    tf_guess,
    free_vars,
    constraints,
    sym_period_targ=1 / 4,
    JCd=None,
    step_size=-1e-2 * 8,
    num_fam_members=40,
    param_continue="jc",
    line_search=True,
)
targeted_po_fam.extend(targeted_po_fam_updated)
for keys in targeted_po_char_updated.keys():
    targeted_po_char[keys].extend(targeted_po_char_updated[keys])

"""
Plot family
"""
# if targeted_po_char != None:
colourby = targeted_po_char['jc']
colourmap='plasma'
cb_label = 'JC'
title = 'EM_L2_Vertical_family'
data_trace = []
# Add L2
data_trace.append(go.Scatter3d(x=[li[0]], y=[0], z=[0], marker=dict(
            color='red',
            size=2)))
# Add Earth
data_trace.append(go.Scatter3d(x=[-mu], y=[0], z=[0], marker=dict(
            color='blue',
            size=10)))
data_trace.append(go.Scatter3d(x=[1-mu], y=[0], z=[0], marker=dict(
            color='grey',
            size=7)))


plot_orbits(mu,targeted_po_fam,colourby, cb_label, title=title,data_trace=data_trace,save=False)
