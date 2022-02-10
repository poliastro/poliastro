"""
Created on Tue Feb  8 22:13:15 2022

@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University
        dhruvj9922@gmail.com

This is an exmaple file to test the CR3BP functions
"""

import numpy as np
import matplotlib.pyplot as plt
from cr3bp_master import prop_cr3bp
from cr3bp_char_quant import bodies_char
from cr3bp_lib_JC_calc import lib_pt_loc, JC

mu, dist_e_m, tstar = bodies_char('Earth','Moon')
print('Earth-Moon mu:',mu)
print('Earth-Moon l*:',dist_e_m,'km')
print('Earth-Moon t*:',tstar/86400,'days')

lib_loc = lib_pt_loc(mu)
li = lib_loc[:,:] # 0 for L1 and  1 for L2....
print('Earth-Moon Li:',li)

ic = [1.05903, -0.067492, -0.103524, -0.170109, 0.0960234, -0.135279]
JC0 = JC(mu,ic[0:3],ic[3:6])
print('Jacobi constant:',JC0)

tf = 10
results = prop_cr3bp(mu, ic, tf,stm_bool = 0,xcross_cond=1)

pltnum = 1
plt.figure(pltnum)
ax = plt.axes(projection='3d')
ax.set_title('CR3BP EM, trajectory, T = 10[nd], tol = 1e-12')
ax.plot3D(results['states'][:,0],results['states'][:,1],results['states'][:,2])
ax.scatter(results['states'][0,0],results['states'][0,1],results['states'][0,2],color='black',label='t=0')
ax.scatter(li[:,0],li[:,1],li[:,2],color='red',label='Li')
ax.scatter(-mu,0,0,color='blue',label='Earth')
ax.scatter(1-mu,0,0,color='grey',label='Moon')
ax.set_box_aspect([ub - lb for lb, ub in (getattr(ax, f'get_{a}lim')() for a in 'xyz')])
ax.set_ylabel("y [nd]")
ax.set_xlabel("x [nd]")
ax.set_zlabel("z [nd]")
plt.legend()
pltnum = pltnum +1 