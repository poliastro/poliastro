"""
@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University
        dhruvj9922@gmail.com

Obj: To compute family of L4 Axial Orbit 
    Single Shooter Variabl Time Setup
    1. Continue in 'vy' + Periodicity (sym_period_targ=1)
    2. Continue in 'jc' + Periodicity (sym_period_targ=1)
    3. Continue in 't' + Periodicity (sym_period_targ=1)

Initial Condition obtained from: 
D. Grebow, "Generating Periodic Orbits in the Circular Restricted Three-Body Problem with Applications to Lunar South Pole Coverage," M.S., May 2006.
"""

import numpy as np
import plotly.graph_objs as go

from cr3bp_char_quant import bodies_char
from cr3bp_lib_JC_calc import lib_pt_loc
from cr3bp_PO_targeter import po_single_shooter_cr3bp
from cr3bp_fam_continuation import npc_po_fam_cr3bp
from cr3bp_plot_orbits import plot_orbits


mu, dist_e_m, tstar = bodies_char("Earth", "Moon")
lib_loc = lib_pt_loc(mu)
li = lib_loc[3,:] # 0 for L1 and  1 for L2

# From D. Grebow
ig = np.array([0.9522, -0.0146, 0.1, -0.1089, 0.0426, 0.4599])
tf_guess = 2.1684

orbit_results = []

free_vars = ['x','y','z','vx','vy','vz','t']
constraints = ['x','y','z','vx','vy'] #vz not constrained as JC implicitly constaints 5 of the 6 states

# Target Axial orbit using Single Shooter Variable Time setup 
#        Exploits Periodicity(sym_perioid_targ set to 1)
#        Continue in 'vy' using Natural Paramter Continuaton 
targeted_po_fam, targeted_po_char = npc_po_fam_cr3bp(mu, po_single_shooter_cr3bp, ig, tf_guess, 
                                      free_vars, constraints, sym_period_targ=1, JCd = None, 
                                      step_size = 1e-3, num_fam_members = 3, param_continue='vy', line_search=True)

# Target Axial orbit using Single Shooter Variable Time setup 
#        Exploits Periodicity(sym_perioid_targ set to 1)
#        Continue in 'jc' using Natural Paramter Continuaton 
ig = targeted_po_fam[-1]['states'][0,:]
tf_guess = targeted_po_fam[-1]['t'][-1]
targeted_po_fam_updated, targeted_po_char_updated = npc_po_fam_cr3bp(mu, po_single_shooter_cr3bp, ig, tf_guess, 
                                      free_vars, constraints, sym_period_targ=1, JCd = None, 
                                      step_size = -1e-3*4, num_fam_members = 150, param_continue='jc', line_search=True)
targeted_po_fam.extend(targeted_po_fam_updated)
for keys in targeted_po_char_updated.keys():
    targeted_po_char[keys].extend(targeted_po_char_updated[keys])

# Target Axial orbit using Single Shooter Variable Time setup 
#        Exploits Periodicity(sym_perioid_targ set to 1)
#        Continue in 't' using Natural Paramter Continuaton 
ig = targeted_po_fam[-1]['states'][0,:]
tf_guess = targeted_po_fam[-1]['t'][-1]
targeted_po_fam_updated, targeted_po_char_updated = npc_po_fam_cr3bp(mu, po_single_shooter_cr3bp, ig, tf_guess, 
                                      free_vars, constraints, sym_period_targ=1, JCd = None, 
                                      step_size = 1e-3*8, num_fam_members = 60, param_continue='t', line_search=True)
targeted_po_fam.extend(targeted_po_fam_updated)
for keys in targeted_po_char_updated.keys():
    targeted_po_char[keys].extend(targeted_po_char_updated[keys])
    
"""
Plot family
"""
if targeted_po_char != None:
    colourby = targeted_po_char['jc']
    colourmap='plasma'
    cb_label = 'JC'
    title = 'EM_L4_Northern_Axial_family'
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
    
    
    plot_orbits(mu,targeted_po_fam,colourby, cb_label, title=title,data_trace=data_trace, save=False)
