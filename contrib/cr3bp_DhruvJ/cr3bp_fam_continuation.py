"""
Created on Mon Feb 21 20:37:16 2022

@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University
    dhruvj9922@gmail.com

Objective: This file contains functions required to continue along a family of Periodic Orbits in the Circular Restricted
    Three Body Problem (CR3BP) model

    Features:
        1. Natural Parameter Continuation (states, time, JC)
        2. Pseduo-Arc Length Continuation

References
____________
This work heavily relies on the work done by the various past and current members of the Multi-Body Dynamics Research Group and Prof. Kathleen C. Howell
These are some of the referneces that provide a comprehensive brackground and have been the foundation for the work:
1. E. Zimovan, "Characteristics and Design Strategies for Near Rectilinear Halo Orbits Within the Earth-Moon System," M.S., August 2017
2. E. Zimovan Spreen, "Trajectory Design and Targeting for Applications to the Exploration Program in Cislunar Space," Ph.D., May 2021
3. V. Szebehely, "Theory of Orbits: The Restricted Problem of Three Bodies", 1967
4. W. Koon, M. Lo, J. Marsden, S. Ross, "Dynamical Systems, The Three-Body Problem, and Space Mission Design", 2006
"""
import copy
import numpy as np
import scipy as sci
from cr3bp_lib_JC_calc import JC
from cr3bp_PO_targeter import map_vars_index_cr3bp

# Add the function as class to inherit from a an orbit class in future to reduce function args
def npc_po_fam_cr3bp(mu, shooter_func, initial_guess, tf_guess, free_vars, constraints, sym_period_targ=1/2, JCd = None, conv_tol=1e-12, int_tol=1e-12, Nmax=10, step_size = 1e-4, num_fam_members = 1, param_continue='x', line_search=False):
    # Can give untargeted member, for PALC need targeted in put
    # line search
    # jc explicity added if passed
    #iterflag = True, False, None

    # Check if paramter to be continued in is defined as a free variable
    if param_continue in free_vars or param_continue == 'jc' :
        free_vars = copy.copy(free_vars)
        constraints = copy.copy(constraints)
        
        # Remove paramter to be continued in from free variable, if not 'jc'
        if param_continue != 'jc':
            free_vars.remove(param_continue)
        
        # To add 'jc' to constraints if family to be continued in JC without explicitly defining it
        if param_continue == 'jc' and 'jc' not in constraints:
            constraints.append('jc')

        param_conti_index = map_vars_index_cr3bp([param_continue])[0]
        
    else:
        print('Paramter that is to be continued in is not defined as a free variable or constraint. Make sure to that the parameter to be continued in can be varied and included in free_vars/constraints')
        return None, None
    
    # Assign a value of JCd if continuing in JC but JCd is not given
    if 'jc' in constraints and JCd == None:
        JCd = JC(mu,initial_guess[0:3],initial_guess[3:6])
        
    targeted_po_fam = []
    targeted_po_char = {'ic':[],'tf':[],'jc':[],'eigenvalues':[],'eigenvectors:':[],'monodromy':[]}
    
    iterflag = False
    count_fam_member = 0
    step_size0 = step_size
    
    while count_fam_member < num_fam_members and iterflag == False:    
        
        results, iterflag = shooter_func(mu, initial_guess, tf_guess, free_vars, constraints, sym_period_targ=sym_period_targ, JCd = JCd, conv_tol=conv_tol, int_tol=int_tol, Nmax=Nmax) 
        
        if iterflag == True:
            # Use Line Search: Update Step size and recompute
            if line_search == True:
                step_size = step_size*0.8
                print('Line search is used to update step size to:',step_size,'\n')
                if param_conti_index < 6:
                    initial_guess[param_conti_index] -= step_size 
                elif param_continue == 't':
                    tf_guess -= step_size
                elif param_continue == 'jc':
                    JCd -= step_size 
                
                if abs(step_size) < abs(step_size0*0.1):
                    print('Updated step size is too small compared to given step size. Rerun with smaller step size')
                else:
                    iterflag = False
                
        elif iterflag == None:
            print('Recheck targeter setup')
            break
        else:
            print('# PO family member = ',count_fam_member+1,'\n')    
            targeted_po_fam.append(results)
            tf_guess = results['t'][-1]
            initial_guess = copy.copy(results['states'][0,:]) # To not update save data as values are passed as object reference
            if param_conti_index < 6:
                initial_guess[param_conti_index] += step_size 
            elif param_continue == 't':
                tf_guess += step_size
            elif param_continue == 'jc':
                JCd += step_size 
            
            # Save key characterisitcs
            targeted_po_char['ic'].append(copy.copy(results['states'][:,0]))
            targeted_po_char['tf'].append(copy.copy(results['t'][-1]))
            targeted_po_char['jc'].append(copy.copy(JC(mu,results['states'][0:3,0],results['states'][3:6,0])))
            targeted_po_char['monodromy'].append(results['stm'][:,:,-1])
            eigenvals, eigenvects = np.linalg.eig(results['stm'][:,:,-1])
            targeted_po_char['eigenvalues'].append(eigenvals)
            targeted_po_char['eigenvectors:'].append(eigenvects)

            
            count_fam_member += 1
            
    return targeted_po_fam, targeted_po_char

def palc_po_fam_cr3bp(mu, shooter_func, targeted_orbit, free_vars, constraints, sym_period_targ=1/2, conv_tol=1e-12, int_tol=1e-12, Nmax=10, step_size = 1e-4, num_fam_members = 1, line_search=False):   

    print('\nAssumes the details of the orbit passed are that of a targeted Periodic Orbit\n')

    if 'jc' in constraints:
        print('JC cannot be constrained when using PALC')
        return None, None
    if len(free_vars) != len(constraints)+1:
        print('Recheck Free variable and constraint setup as Null space needs to be exactly one')
        return None, None

    # Setup PALC arguments to target PALC based orbits
    palc_args = {}    
    palc_args['delta_s'] = step_size
    
    # Compute Null Space
    xfree_null_prev = sci.linalg.null_space(targeted_orbit['DF'])
    if np.size(xfree_null_prev, 1) != 1:
        print('Null space is not one, nullity is',np.size(xfree_null_prev, 1))
        break
    
    xfree_null_prev = xfree_null_prev.flatten()
    null_vecs_dot = np.dot(xfree_null_prev, null_vec)
    null_vec = xfree_null_prev*np.sign(null_vecs_dot)

    palc_args['xfree_prev'] = results_stm['xfree']
    palc_args['delta_X*_prev'] = null_vec#/xfree_null_prev[0,0]
    
    return targeted_po_fam, targeted_po_char