"""
Created on Mon Feb 21 20:37:16 2022

@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University
    dhruvj9922@gmail.com

Objective: This file contains functions required to target a Periodic Orbit in the Circular Restricted
    Three Body Problem (CR3BP) model

    Features:
        1. Setup multiple nodes to use a Multiple Shooter Targeter
        2. (will be added)Periodic Orbit Multiple Shooter Targeter (Acts as a single shooter if n_node == 1)
        3. (will be added)Newton-Raphson Solver

References
____________
This work heavily relies on the work done by the various past and current members of the Multi-Body Dynamics Research Group and Prof. Kathleen C. Howell
These are some of the referneces that provide a comprehensive brackground and have been the foundation for the work:
1. E. Zimovan, "Characteristics and Design Strategies for Near Rectilinear Halo Orbits Within the Earth-Moon System," M.S., August 2017
2. E. Zimovan Spreen, "Trajectory Design and Targeting for Applications to the Exploration Program in Cislunar Space," Ph.D., May 2021
3. V. Szebehely, "Theory of Orbits: The Restricted Problem of Three Bodies", 1967
4. W. Koon, M. Lo, J. Marsden, S. Ross, "Dynamical Systems, The Three-Body Problem, and Space Mission Design", 2006
"""
import numpy as np
from cr3bp_lib_JC_calc import JC
from cr3bp_master import prop_cr3bp, ui_partials_acc_cr3bp

def po_single_shooter_cr3bp(mu, initial_guess, tf_guess, free_vars, constraints, sym_period_targ=1/2, palc_args=None, conv_tol=1e-12, int_tol=1e-12, Nmax=50):
    # sym_period_targ = 1/2 -> symmetry after half period, 1/4 -> after quarter period, 1-> periodicity
    # lyap, halo, axial, vertical write free_vars and constraints for each and sym_period_targ
    # tfguess = period
    # xd
    #palc
    # note care about palc and JC rn
    # test periodicity
    # To do: palc, JC constraint, periodicity and xd
    
    if 'jc' in free_vars:
        print('Jacobi Constant cannot be a free variable')
        return 0
    if 't' in constraints:
        print('Time should not be a constraint, instead make the tf_guess to be desired time and not add \'t\' as a free variable')
    
    if sym_period_targ not in [1/4,1/2,1]:
        print('Not a valid fraction of period to target')
        return 0
    
    free_vars_index, constraints_index = map_vars_index_cr3bp(free_vars, constraints)
        
    # Use events function to compute ycrossing if symmetry is to be used for targeter
    if ('t' not in free_vars or sym_period_targ !=1):        
        results_stm = prop_cr3bp(mu, initial_guess, tf_guess, tol=int_tol, xcross_cond=1)
        if len(results_stm['tevents'][0][:]) > 1: # condition so that it works for corrected solutions
            tf = results_stm['tevents'][0][1]
        else:    # Assumes that IC doesnt count as first pass, the next crossing that is needed for targeting is the yevent
            tf = results_stm['tevents'][0][0]
    else:# If need to used time fixed targeter or periodicity targeter 
        tf = tf_guess*sym_period_targ 

    results_stm = prop_cr3bp(mu, initial_guess, tf, stm_bool=1, tol=int_tol, xcross_cond=0)        
    ic = initial_guess
    
    xfree = np.zeros(len(free_vars_index))
    xconstraint = np.zeros(len(constraints_index))
    xdesired = np.zeros(len(constraints_index))
    DF = np.zeros((len(xfree),len(xconstraint)))
    
    print(constraints_index, free_vars_index, tf)
    
    stm_col_index = [free_vars_index[i] for i in range(len(xfree)) if free_vars_index[i] < 6]
    stm_row_index = [constraints_index[i] for i in range(len(xconstraint)) if constraints_index[i] < 6]
    stm_col_len = len(stm_col_index)
    stm_row_len = len(stm_row_index)
    
    if 't' in free_vars:
        xfree[:-1] = ic[stm_col_index]
        xfree[-1] = tf
    else:
        xfree = ic[stm_col_index]
    
    if 'jc' in constraints:
        xconstraint[:-1] = results_stm['states'][-1,stm_row_index]
        xconstraint[-1] = JC(mu,ic[0:3],ic[3:6])
    else:
        xconstraint = results_stm['states'][-1,stm_row_index]

    count = 0
    
    iterflag = False
    
    # Targeter loop: Use L2-norm of constraint vector and #iterations as stopping condition
    while np.linalg.norm(xconstraint-xdesired) > conv_tol and count <= Nmax:
    
        # Update STM after each iteration to converge faster
        states_final = results_stm["states"][-1,:]
        Ux, Uy, Uz, ax, ay, az = ui_partials_acc_cr3bp(mu, states_final)
        Dot_states = np.array([states_final[3], states_final[4], states_final[5], ax, ay, az])
        
        stm_temp = results_stm["stm"][stm_row_index, :, -1]
        DF[:stm_row_len,:stm_col_len] = stm_temp[:, stm_col_index]  
        
        if 't' in free_vars:
            DF[:,-1] = Dot_states[stm_row_index] 
            
        FX = xconstraint-xdesired
        xfree = newton_raphson_update(xfree, FX, DF)
        
        if 't' in free_vars:
            ic[stm_col_index] = xfree[:-1]
            tf = xfree[-1]
        else:
            ic[stm_col_index] = xfree

        results_stm = prop_cr3bp(mu, ic, tf, stm_bool=1, tol=int_tol, xcross_cond=0)        
        
        if 'jc' in constraints:
            xconstraint[:-1] = results_stm['states'][-1,stm_row_index]
            xconstraint[-1] = JC(mu,ic[0:3],ic[3:6])
        else:
            xconstraint = results_stm['states'][-1,stm_row_index]
        
        print(count, np.linalg.norm(xconstraint-xdesired), np.argmax(np.abs(xconstraint-xdesired)),np.max(np.abs(xconstraint-xdesired)))
        
        count = count + 1 
    
    if count > Nmax:
        iterflag = True
        print('STOPPPPPPPPPPPPPPPPPPPP')
        
    results_stm = prop_cr3bp(mu, ic, tf*1/sym_period_targ, stm_bool=1, tol=int_tol, xcross_cond=0)      
    
    return results_stm, iterflag

def newton_raphson_update(xfree,FX,DG):
    """
    Multi-dimensional Newton-Raphson Method to update inital guess
    Dhruv Jain, Feb 26 2022

    Parameters
    ----------
    xfree : numpy ndarray, float
        Free variables
    FX : numpy ndarray, float
        Constraint vector
    DG : numpy ndarray, float
        Jacobian Matrix
        
    Returns
    -------
    xfree : numpy ndarray, flaot
        Updated Free variables
    """
    
    if len(xfree) == len(FX):
        DG_inv = np.linalg.inv(DG)
    else: # Not a square matrix, use Pseudo-inverse
        DG_inv = np.linalg.pinv(DG)
    
    delta_x_guess = np.matmul(DG_inv, FX)    
    xfree = xfree - delta_x_guess
    
    return xfree

def multi_shooter_nodes_setup_cr3bp(mu, ic, tf, n_node, node_place_opt="time", int_tol=1e-12):
    """
    Calculate states of n_nodes and time between each node
    Patch point placement strategy:
        1) Computes the n_nodes placed after NEARLY equal time intervals
        2) Computes the n_nodes placed after NEARLY equal time history INDEX
            (This might be better as more integration steps are taken in sensitive
             regions and placing the nodes by using index will place more points
             in the sensitive regions)

    First Node = I.C. of trajectory
    If symmetry is to be leverage to target a P.O.: last Node propagated by time 't' reaches the vicinity of the desired final state

    Dhruv Jain, Feb 18 2022

    Parameters
    ----------
    mu :  float, M2/(M1+M2)
        M1 and M2 are mass of Primary Bodies and M2<M1
    ic : numpy ndarray (6x1), {Can handle all 42 states for CR3BP+STM integration}
        States are defined about the barycenter of the two primaries, P1 and P2
        Initial condition: 6 states to compute a trajectory];
        [0:x0, 1:y0, 2:z0, 3:vx0, 4:vy0, 5:vz0] [non-dimensional] [nd]
    tf : float
        Integration time [nd]
        Can be negative or positive, negative => Integration in backwards time
    n_node : int, Number of nodes
        IC is node 1 and nth node is node when propagated by tn value should reach the vicinity of the desired final state
    node_place_opt : string, optional
        'time': Computes the n_nodes placed after NEARLY equal time intervals
        'index': Computes the n_nodes placed after NEARLY equal time history INDEX
    int_tol : float, optional
        Absolute = Relative Integration Tolerance
        The default is 1e-12.

    Returns
    -------
    ic_node : list, ndarray, float64
        Stores the IC of the n_node, where  ic_node[0] = ic0
    t_node : list, float
        Stores the time from one node to the next, t_node[-1]: time to reach the desired state (somekind of corrsing)
    """

    if node_place_opt != "time" and node_place_opt != "index":
        print(
            "Incorrect node placement option passed. Allowable options: time and index"
        )
        return 0

    results_stm = prop_cr3bp(
        mu, ic, tf, stm_bool=0, xcross_cond=0
    )  # Propagate I.C. till tf
    ic_node = []
    t_node = []

    ic_node.append(ic)

    if node_place_opt == "time":
        # Time of each segment
        ti = np.linspace(0, tf, n_node + 1)
        ti = ti[
            1:
        ]  # As ti is the time from node_i to next node, omit the first time, i.e. 0

        for i in range(1, n_node):
            index = np.argmin(
                abs(results_stm["t"] - ti[i - 1])
            )  # compute index when index from time history when time is nearly
            ic_node.append(results_stm["states"][index, :])
            t_node.append(
                results_stm["t"][index] - sum(t_node)
            )  # t_node runs one node behind as it is the time to next node, subtract sum as ti+t2+t3 = tf

        t_node.append(
            tf - sum(t_node)
        )  # Time from final node to vicinity of desired state

    elif node_place_opt == "index":
        num_index = len(results_stm["t"])
        indices = np.linsapce(
            0, num_index, n_node, dtype="int"
        )  # Linearly spaced indices

        for i in range(1, n_node):
            ic_node.append(results_stm["states"][indices[i], :])
            t_node.append(
                results_stm["t"][indices[i]] - sum(t_node)
            )  # t_node runs one node behind as it is the time to next node, subtract sum as ti+t2+t3 = tf

        t_node.append(
            tf - sum(t_node)
        )  # Time from final node to vicinity of desired state

    return ic_node, t_node

def map_vars_index_cr3bp(free_vars, constraints):
    """
    Map variables defined in free variables and constraints to numerical index
    Dhruv Jain, Feb 25 2022
    
    Parameters
    ----------
    free_vars : list of strings
        Contains strings expressing the free variables
    constraints : list of strings
        Contains strings expressing the constraints

    Returns
    -------
    free_vars_index : list of int
        Integer index code to represent free variables
    constraints_index : list of int
        Integer index code to represnet constraints
    """    
    
    variable_dict = {'x': 0, 'y': 1, 'z': 2, 'vx': 3, 'vy': 4, 'vz': 5, 'jc': 6, 't':7}
    
    free_vars_index = []
    constraints_index = []
    
    for i in range(len(free_vars)):
        # Check and handle KeyError
        try: 
            free_vars_index.append(variable_dict[free_vars[i]])
        except KeyError:
            print(free_vars[i],'is not a valid free variable')
            free_vars_index = 0
            break
            
    for j in range(len(constraints)):
        # Check and handle KeyError
        try:
            constraints_index.append(variable_dict[constraints[j]])
        except KeyError:
            print(constraints[j],'is not a valid constraint')
            constraints_index = 0
            break        
    
    # Sort the indices to handle variables passed in any order
    free_vars_index = sorted(free_vars_index)
    constraints_index = sorted(constraints_index)
    
    return free_vars_index, constraints_index