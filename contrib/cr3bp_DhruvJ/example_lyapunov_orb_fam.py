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
from cr3bp_char_quant import sys_chars
from cr3bp_lib_calc import lib_pt_loc
from cr3bp_lyap_ig_generator import ig_lyap_orb_collinear_li_cr3bp
from cr3bp_po_fam_continuation import periodic_orbit_fam_continuation
from cr3bp_po_plot_orbits import plot_orbits

sys_p1p2 = sys_chars("Earth", "Moon")
mu = sys_p1p2.mu
lib_loc = lib_pt_loc(sys_p1p2)
li = lib_loc[1, :]  # 0 for L1 and  1 for L2

# Compute Initial Guess of L1 Lyapunov Orbit using linar approximation
initial_guess = ig_lyap_orb_collinear_li_cr3bp(sys_p1p2, pert_x=0.01)

tf_guess = 5  # Give a random large enough tf as tf_guess, events function with XZ plane symmetry will correct it

# Lyapunov is a planar orbit - [x, y, vx, vy]
free_vars = ["x", "vy", "t"]  # All the possible free variables
constraints = [
    "y",
    "vx",
]  # All constraints, 'jc' can be added explicity or will be added in npc_fam_fam_cr3bp if to be continued in 'jc'

orb_fam_obj = periodic_orbit_fam_continuation(
    sys_p1p2, initial_guess, tf=tf_guess
)

# Target Lyapunov orbit using Single Shooter Variable Time setup
#         Exploits XZ plane symmetry(sym_perioid_targ set to 1/2)
#         Continue in 'x' using Natural Paramter Continuaton to compute 10 family members
orb_fam_obj.npc_po_fam(
    free_vars,
    constraints,
    sym_period_targ=1 / 2,
    Nmax=20,
    step_size=1e-4,
    num_fam_members=10,
    param_continue="x",
    line_search=False,
)


# Target Lyapunov orbit using Single Shooter Variable Time setup
#        Exploits Periodcity(sym_perioid_targ set to 1)
#        Continue in 'JC' using Natural Paramter Continuaton to compute 10 family members
# Note: JC coonstraint is added without being explicity defined as a constraint
orb_fam_obj.JCd = orb_fam_obj.JC() - 1e-3
orb_fam_obj.npc_po_fam(
    free_vars,
    constraints,
    sym_period_targ=1,
    Nmax=20,
    step_size=-1e-3 * 5,
    num_fam_members=10,
    param_continue="jc",
    line_search=False,
)


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
orb_fam_obj.JCd -= 1e-5
orb_fam_obj.npc_po_fam(
    free_vars,
    constraints,
    sym_period_targ=1,
    Nmax=10,
    step_size=-1e-3 * 5,
    num_fam_members=20,
    param_continue="jc",
    line_search=True,
)

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

orb_fam_obj.JCd = orb_fam_obj.JC()
orb_fam_obj.npc_po_fam(
    free_vars,
    constraints,
    sym_period_targ=1 / 2,
    Nmax=10,
    step_size=-1e-3 * 5,
    num_fam_members=10,
    param_continue="jc",
    line_search=True,
)

"""
Plot family
"""
colourby = orb_fam_obj.targeted_po_char["jc"]
colourmap = "plasma"
cb_label = "JC"
title = "EM_L1_Lyapunov_family"
data_trace = []
data_trace.append(
    go.Scatter3d(x=[1 - mu], y=[0], z=[0], marker=dict(color="grey", size=7))
)

plot_orbits(
    mu,
    orb_fam_obj.targeted_po_fam,
    colourby,
    cb_label,
    title=title,
    save=False,
)
