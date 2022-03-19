"""
@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University
        dhruvj9922@gmail.com

Obj: To compute family of L1 Lyapunov Orbit
    Test linear initial guess for Lyapunov Orbit
    Single Shooter Variabl Time Setup
    1. Continue in 'x' + XZ plane symmetry use => targets Period/2 states
    2. Continue in 'jc' + Periodicity targeter => targets Period states
    3. Continue in 'jc' + Periodicity targeter with line_search=> targets Period states
    4. Continue in 'jc' + XZ plane symmetry targeter with line_search=> targets Period states

Note: If the step size is too big and targeting periodicity with FX = ['y','vx']then the may converge
to states at next XZ plane corssing instead of targeting states after 1 period
"""
import plotly.graph_objs as go
from cr3bp_char_quant import bodies_char
from cr3bp_fam_continuation import npc_po_fam_cr3bp
from cr3bp_initial_guess_generator import ig_lyap_orb_collinear_li_cr3bp
from cr3bp_lib_JC_calc import JC, lib_pt_loc
from cr3bp_plot_orbits import plot_orbits
from cr3bp_PO_targeter import po_single_shooter_cr3bp

mu, dist_e_m, tstar = bodies_char("Earth", "Moon")
lib_loc = lib_pt_loc(mu)
li = lib_loc[0, :]  # 0 for L1 and  1 for L2

# Compute Initial Guess of L1 Lyapunov Orbit using linar approximation
initial_guess = ig_lyap_orb_collinear_li_cr3bp(mu, pert_x=0.01)

tf_guess = 5  # Give a random large enough tf as tf_guess, events function with XZ plane symmetry will correct it

# Lyapunov is a planar orbit - [x, y, vx, vy]
free_vars = ["x", "vy", "t"]  # All the possible free variables
constraints = [
    "y",
    "vx",
]  # All constraints, 'jc' can be added explicity or will be added in npc_fam_fam_cr3bp if to be continued in 'jc'

# Target Lyapunov orbit using Single Shooter Variable Time setup
#        Exploits XZ plane symmetry(sym_perioid_targ set to 1/2)
#        Continue in 'x' using Natural Paramter Continuaton to compute 10 family members
targeted_po_fam, targeted_po_char = npc_po_fam_cr3bp(
    mu,
    po_single_shooter_cr3bp,
    initial_guess,
    tf_guess,
    free_vars,
    constraints,
    sym_period_targ=1 / 2,
    JCd=None,
    step_size=1e-4,
    num_fam_members=10,
    param_continue="x",
    line_search=False,
)

# Target Lyapunov orbit using Single Shooter Variable Time setup
#        Exploits Periodcity(sym_perioid_targ set to 1)
#        Continue in 'JC' using Natural Paramter Continuaton to compute 10 family members
# Note: JC coonstraint is added without being explicity defined as a constraint
initial_guess = targeted_po_fam[-1]["states"][0, :]
tf_guess = targeted_po_fam[-1]["t"][-1]
JCd = JC(mu, initial_guess[0:3], initial_guess[3:6]) - 1e-3
targeted_po_fam_jc, targeted_po_char_jc = npc_po_fam_cr3bp(
    mu,
    po_single_shooter_cr3bp,
    initial_guess,
    tf_guess,
    free_vars,
    constraints,
    sym_period_targ=1,
    JCd=JCd,
    step_size=-1e-3 * 5,
    num_fam_members=10,
    param_continue="jc",
    Nmax=10,
    line_search=False,
)
# Add the results of 'jc' continuation to 'x' continuation family members
targeted_po_fam.extend(targeted_po_fam_jc)
for keys in targeted_po_char_jc.keys():
    targeted_po_char[keys].extend(targeted_po_char_jc[keys])

"""
Constraint updated to show that periodicity can be targeted with other states + line search
"""
constraints = [
    "y",
    "vy",
]  # All constraints, 'jc' can be added explicity or will be added in npc_fam_fam_cr3bp if to be continued in 'jc'


# Target Lyapunov orbit using Single Shooter Variable Time setup
#        Exploits Periodcity(sym_perioid_targ set to 1)
#        Continue in 'JC' using Natural Paramter Continuaton to compute 30 family members
#        Line search is used to update step size if unable to converge
# Note: JC coonstraint is added without being explicity defined as a constraint
initial_guess = targeted_po_fam[-1]["states"][0, :]
tf_guess = targeted_po_fam[-1]["t"][-1]
JCd = JC(mu, initial_guess[0:3], initial_guess[3:6]) - 1e-3
targeted_po_fam_jc, targeted_po_char_jc = npc_po_fam_cr3bp(
    mu,
    po_single_shooter_cr3bp,
    initial_guess,
    tf_guess,
    free_vars,
    constraints,
    sym_period_targ=1,
    JCd=JCd,
    step_size=-1e-3 * 5,
    num_fam_members=20,
    param_continue="jc",
    Nmax=10,
    line_search=True,
)
# Add the results of 'jc' continuation to 'x' continuation family members
targeted_po_fam.extend(targeted_po_fam_jc)
for keys in targeted_po_char_jc.keys():
    targeted_po_char[keys].extend(targeted_po_char_jc[keys])

"""
Constraint updated to use XZ plane symmetry
"""
constraints = [
    "y",
    "vx",
]  # All constraints, 'jc' can be added explicity or will be added in npc_fam_fam_cr3bp if to be continued in 'jc'

# Target Lyapunov orbit using Single Shooter Variable Time setup
#        Exploits XZ plane crossing (sym_perioid_targ set to 1/2)
#        Continue in 'JC' using Natural Paramter Continuaton to compute 100 family members
#        Line search is used to update step size if unable to converge
# Note: JC coonstraint is added without being explicity defined as a constraint
initial_guess = targeted_po_fam[-1]["states"][0, :]
tf_guess = targeted_po_fam[-1]["t"][-1]
JCd = JC(mu, initial_guess[0:3], initial_guess[3:6])
targeted_po_fam_jc, targeted_po_char_jc = npc_po_fam_cr3bp(
    mu,
    po_single_shooter_cr3bp,
    initial_guess,
    tf_guess,
    free_vars,
    constraints,
    sym_period_targ=1 / 2,
    JCd=JCd,
    step_size=-1e-3 * 5,
    num_fam_members=65,
    param_continue="jc",
    Nmax=10,
    line_search=True,
)
# Add the results of 'jc' continuation to 'x' continuation family members
targeted_po_fam.extend(targeted_po_fam_jc)
for keys in targeted_po_char_jc.keys():
    targeted_po_char[keys].extend(targeted_po_char_jc[keys])

"""
Plot family
"""
colourby = targeted_po_char["jc"]
colourmap = "plasma"
cb_label = "JC"
title = "EM_L1_Lyapunov_family"
data_trace = []
data_trace.append(
    go.Scatter3d(x=[1 - mu], y=[0], z=[0], marker=dict(color="grey", size=7))
)

plot_orbits(mu, targeted_po_fam, colourby, cb_label, title=title, save=False)
