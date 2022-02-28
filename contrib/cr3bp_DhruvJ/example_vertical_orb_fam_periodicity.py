"""
Created on Sat Feb 26 23:19:40 2022

@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University
        dhruvj9922@gmail.com

NPC of Vertical family in vy and JC
"""
import copy
import numpy as np
import matplotlib.pyplot as plt
from cr3bp_char_quant import bodies_char
from cr3bp_lib_JC_calc import lib_pt_loc, JC
from cr3bp_master import prop_cr3bp
from cr3bp_initial_guess_generator import ig_lyap_orb_collinear_li_cr3bp
from cr3bp_PO_targeter import po_single_shooter_cr3bp

mu, dist_e_m, tstar = bodies_char("Earth", "Moon")
lib_loc = lib_pt_loc(mu)
li = lib_loc[0,:] # 0 for L1 and  1 for L2

orbit_results = []

# ig= ig_lyap_orb_collinear_li_cr3bp(mu, pert_x=0.01)

# tf_guess = 5

# free_vars = ['vy','t']
# constraints = ['y','vx']

# results, iterflag = po_single_shooter_cr3bp(mu, ig, tf_guess, free_vars, constraints)
# orbit_results.append(results)
# tf_guess = results['t'][-1]
# ig = copy.copy(results['states'][0,:])
# ig[0] += 0.0008 
    
    
# free_vars = ['vx','vy','t']
# constraints = ['x','y','vx']


ig = np.array([1.0842, 0, 0, 0, -0.5417,0.8415])
tf_guess = 6.1305

free_vars = ['x','vz','t']
constraints = ['x','vy','vz']
    
for i in range(2):
    results, iterflag = po_single_shooter_cr3bp(mu, ig, tf_guess, free_vars, constraints,sym_period_targ=1)
    orbit_results.append(results)
    tf_guess = results['t'][-1]
    print('JC:',JC(mu,results['states'][0,0:3],results['states'][0,3:6]),tf_guess)
    ig = copy.copy(results['states'][0,:])
    print(ig)
    ig[4] += -0.0001

free_vars = ['x','vy','vz','t']
constraints = ['x','vy','vz','jc']
JCd = round(JC(mu,results['states'][0,0:3],results['states'][0,3:6]),4)

for i in range(15):
    results, iterflag = po_single_shooter_cr3bp(mu, ig, tf_guess, free_vars, constraints,JCd=JCd,sym_period_targ=1)
    orbit_results.append(results)
    tf_guess = results['t'][-1]
    print('JC:',JC(mu,results['states'][0,0:3],results['states'][0,3:6]),tf_guess)
    JCd += -4.5e-2

    
plt.figure(1)
ax = plt.axes(projection='3d')
ax.set_title('EM, L1 Halo Orbit Family, tol = 1e-12')
for i in range(len(orbit_results)):
    ax.plot3D(orbit_results[i]["states"][:,0], orbit_results[i]["states"][:,1],orbit_results[i]["states"][:,2])
# plt.plot(li[0],li[1],'ro',label='L1')
# plt.plot(1-mu,0,'b*')
ax.scatter(li[0],li[1],li[2],color='red')
ax.set_box_aspect([ub - lb for lb, ub in (getattr(ax, f'get_{a}lim')() for a in 'xyz')])
ax.set_ylabel("y [nd]")
ax.set_xlabel("x [nd]")
ax.set_zlabel("z [nd]")
plt.show()
