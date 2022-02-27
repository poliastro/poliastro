"""
Created on Sat Feb 26 20:46:54 2022

@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University
        dhruvj9922@gmail.com
        
        
"""
import copy
import matplotlib.pyplot as plt
from cr3bp_char_quant import bodies_char
from cr3bp_lib_JC_calc import lib_pt_loc
from cr3bp_master import prop_cr3bp
from cr3bp_initial_guess_generator import ig_lyap_orb_collinear_li_cr3bp
from cr3bp_PO_targeter import po_single_shooter_cr3bp


mu, dist_e_m, tstar = bodies_char("Earth", "Moon")
lib_loc = lib_pt_loc(mu)
li = lib_loc[0,:] # 0 for L1 and  1 for L2

initial_guess = ig_lyap_orb_collinear_li_cr3bp(mu, pert_x=0.01)

tf_guess = 5

free_vars = ['vy','t']
constraints = ['y','vx']

lyap_results = []

for i in range(30):
    results, iterflag = po_single_shooter_cr3bp(mu, initial_guess, tf_guess, free_vars, constraints)
    lyap_results.append(results)
    tf_guess = results['t'][-1]
    initial_guess = copy.copy(results['states'][0,:])
    initial_guess[0] += 0.0008 


plt.figure(1)
plt.title("L1 Lyapunov Orbit Family : Dhruv Jain")
for i in range(len(lyap_results)):
    plt.plot(lyap_results[i]["states"][:,0], lyap_results[i]["states"][:,1])
plt.plot(li[0],li[1],'ro',label='L1')
# plt.plot(1-mu,0,'b*')
plt.xlabel("x [nd]")
plt.ylabel("y [nd]")
plt.axis('equal')
plt.grid()
plt.show()
