"""
Created on Mon Feb 21 20:37:16 2022

@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University
    dhruvj9922@gmail.com

Objective: This file contains functions required to target a Periodic Orbit in the Circular Restricted
    Three Body Problem (CR3BP) model

    Features:
        1. Single Shooter Targeter: 
            -> Variable time targeter capable of targeting states, JC and periodicity 
            -> Pseudo-Arc Length Compatible
            -> Line Search to update step size         
        2. Newton-Raphson Solver 
        3. Setup multiple nodes to use a Multiple Shooter Targeter

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
import copy
from cr3bp_lib_JC_calc import JC
from cr3bp_master import prop_cr3bp, ui_partials_acc_cr3bp

def po_single_shooter_cr3bp(mu, initial_guess, tf_guess, free_vars, constraints, sym_period_targ=1/2, JCd = None, palc_args=None, conv_tol=1e-12, int_tol=1e-12, Nmax=50):
    """Single shooter targeter for a Periodic Orbit defined in the CR3BP model
    Dhruv Jain, 26 Feb 2022

    Parameters
    ----------
    mu :  float, M2/(M1+M2)
        M1 and M2 are mass of Primary Bodies and M2<M1
    initial_guess : numpy ndarray (6x1)
        States are defined about the barycenter of the two primaries, P1 and P2
        Initial guess: 6 states to target a Periodic Orbti;
        [0:x0, 1:y0, 2:z0, 3:vx0, 4:vy0, 5:vz0] [non-dimensional] [nd]
    tf_guess : float
        Period guess of a periodic orbit[nd], should be positive, not validated for negative time
    free_vars : list of strings
        Describes the parameters set as free variables
        Possible paramters: ['x','y','z','vx','vy','vz','t']
    constraints : list of strings
        Describes the parameters set as constraints
        Possible paramters: ['x','y','z','vx','vy','vz','jc']
    sym_period_targ : float, optional
        Describes the fraction of period to be targeted. The default is 1/2.
        1/4: Usually use for vertical orbits to use XZ plane AND X-axis symmetry to target states after Period/4
        1/2: Usually use for lyapunov, halo, axial to leverage XZ plane OR X-axis symmetry to target states after Period/2
        1: Usually use to target orbits that do not have XZ plane OR X-axis symmetry, target orbit periodicity
    JCd : float, optional
        Desired value of JC to be targeted. The default is None.
    palc_args : dictionary containting information for Pseudo-Arc Length Continuation, optional
        'delta_s': step-size
        'free_var_prev': Free variables of the previously converged orbit
        'delta_X*_prev': Null-Vector of DF of previously converged solution
        The default is None
    conv_tol : float, optional
        Convergence tolerance of the constraint vector -> Acceptable tolerance of L2 norm of the constraint vector . The default is 1e-12.
    int_tol : float, optional
        Absolute = Relative Integration Tolerance
        The default is 1e-12.
    Nmax : int, optional
        Max allowable iteration of targeter while loop. The default is 50.

    Returns
    -------
    results : Targeted Periodic Orbit - Dictionary
        't': time history, [nd]
        'states': state history, [:,i] => 6 states @ t = t_i, [nd]
        'yevents': states at prescribed event, [nd]
        'tevents': time stamp of yevents, [nd]
        'stm': STM element history, [:,:,i] => STM @ t = t_i
    iterflag: Boolean/None
        True: Targeter unable to converge
        None: Targeter Setup is incorrect
        False: Succesfully targeted
    """
    
    iterflag = None
    # Free vairables and Constraints paramters logical error check    
    if 'jc' in free_vars:
        print('Jacobi Constant cannot be a free variable')
        return 0, iterflag
    if 't' in constraints:
        print('Time should not be a constraint, instead make the tf_guess to be desired time and not add \'t\' as a free variable')
        return 0, iterflag
    if sym_period_targ not in [1/4,1/2,1]:
        print('Not a valid fraction of period to target')
        return 0
    
    if palc_args != None:
        print('PALC constraint included') 
    
    # Map paramter strings to index
    free_vars_index = map_vars_index_cr3bp(free_vars)
    constraints_index = map_vars_index_cr3bp(constraints)
        
    
    # Use events function to compute y-crossing if symmetry is to be used for targeter
    if 't' in free_vars and sym_period_targ !=1:        
        results_stm = prop_cr3bp(mu, initial_guess, tf_guess, tol=int_tol, xcross_cond=1)
        if len(results_stm['tevents'][0][:]) > 1: # condition so that it works for corrected solutions
            tf = results_stm['tevents'][0][1]
        else:    # Assumes that IC doesnt count as first pass, the next crossing that is needed for targeting is the yevent
            tf = results_stm['tevents'][0][0]
    else:# If need to used time fixed targeter or periodicity targeter 
        tf = tf_guess*sym_period_targ 
    
    results_stm = prop_cr3bp(mu, initial_guess, tf, stm_bool=1, tol=int_tol, xcross_cond=0)        
    ic = copy.copy(initial_guess)
    
    
    # Initialize key components of targeter
    xfree = np.zeros(len(free_vars_index))
    xconstraint = np.zeros(len(constraints_index))
    xdesired = np.zeros(len(constraints_index))
    DF = np.zeros((len(xconstraint),len(xfree)))
    DG = np.zeros((len(xfree),len(xfree))) # DF for PALC
    
    print('FX:', constraints,'X:', free_vars,'\nTf0:', tf,'JC0:', JC(mu,ic[0:3],ic[3:6]))
    
    stm_col_index = [free_vars_index[i] for i in range(len(xfree)) if free_vars_index[i] < 6]
    stm_row_index = [constraints_index[i] for i in range(len(xconstraint)) if constraints_index[i] < 6]
    stm_col_len = len(stm_col_index)
    stm_row_len = len(stm_row_index)
    
    # Setup Free Variable Vector
    if 't' in free_vars:
        xfree[:-1] = ic[stm_col_index]
        xfree[-1] = tf
    else:
        xfree = ic[stm_col_index]
    
    # Setup Constraint Vector, xconstraint where FX = xconstraint - xdesired
    if 'jc' in constraints:
        xconstraint[:-1] = results_stm['states'][-1,stm_row_index]
        xconstraint[-1] = JC(mu,ic[0:3],ic[3:6])
    else:
        xconstraint = results_stm['states'][-1,stm_row_index]
    
    # Setup Desired vetor, xdesired where FX = xcontraint - xdesired
    if 'jc' in constraints:
        # xdesired[:-1] = 0 #results_stm['states'][0,stm_row_index]  
        xdesired[-1] = JCd
    # Updated desired to be inital state of orbit    
    if sym_period_targ == 1:
        xdesired[:stm_row_len] = results_stm['states'][0,stm_row_index]    
    if sym_period_targ == 1:
        # Create identity like matrix to be subtracted from STM    
        identity_mat = np.eye(6)
        temp = identity_mat[stm_row_index,:]
        identity_temp = temp[:,stm_col_index]
    
    # Create FX        
    if palc_args == None:        
        FX = xconstraint-xdesired
    else:
        FX = np.zeros((len(xconstraint)+1))
        FX[:-1] = xconstraint-xdesired
        FX[-1] = np.matmul(np.transpose(xfree-palc_args['free_var_prev']), palc_args['delta_X*_prev']) - palc_args['delta_s'] 
        
    count = 0
        
    print('Iteration:',count,'|FX|=',np.linalg.norm(FX))
   
    iterflag = False    
    # Targeter loop: Use L2-norm of constraint vector and #iterations as stopping condition
    while np.linalg.norm(FX) > conv_tol and count <= Nmax:    
        # Update STM after each iteration to converge faster
        
        states_final = results_stm["states"][-1,:]
        Ux, Uy, Uz, ax, ay, az = ui_partials_acc_cr3bp(mu, states_final)
        Dot_states = np.array([states_final[3], states_final[4], states_final[5], ax, ay, az])
        
        # Compute d(JC)/d(x)
        Ux_ic, Uy_ic, Uz_ic, _,_,_ = ui_partials_acc_cr3bp(mu, ic)
        dJC_dx = np.array([2*Ux_ic, 2*Uy_ic, 2*Uz_ic, -2*ic[3], -2*ic[4], -2*ic[5]])
        
        #Extract required elements of the Monodromy matrix for the targeter setup
        stm_temp = results_stm["stm"][stm_row_index, :, -1]


        # Setup DF, Jacobian Matrix
        DF[:stm_row_len,:stm_col_len] = stm_temp[:, stm_col_index]          
        if 'jc' in constraints:
            DF[-1,:stm_col_len] = dJC_dx[stm_col_index]        
        if 't' in free_vars:
            DF[:stm_row_len,-1] = Dot_states[stm_row_index] 
        if sym_period_targ == 1:
            DF[:stm_row_len,:stm_col_len] = DF[:stm_row_len,:stm_col_len] - identity_temp           
            
        # Update Free variable vector, include PALC constraint if PALC is being used        
        if palc_args == None:        
            xfree = newton_raphson_update(xfree, FX, DF)
        else:
            DG[:-1,:] = DF
            DG[-1,:] = palc_args['delta_X*_prev'] 
            xfree = newton_raphson_update(xfree, FX, DG)
        
        # Update Initial Condition and Time
        if 't' in free_vars:
            ic[stm_col_index] = xfree[:-1]
            tf = xfree[-1]
        else:
            ic[stm_col_index] = xfree

        results_stm = prop_cr3bp(mu, ic, tf, stm_bool=1, tol=int_tol, xcross_cond=0)        
        
        # Update xconstraint
        if 'jc' in constraints:
            xconstraint[:-1] = results_stm['states'][-1,stm_row_index]
            xconstraint[-1] = JC(mu,ic[0:3],ic[3:6])
        else:
            xconstraint = results_stm['states'][-1,stm_row_index]

        # Update xdesired            
        if 'jc' in constraints:
            xdesired[-1] = JCd            
            # Updated desired to be inital state of orbit    
        if sym_period_targ == 1:
            xdesired[:stm_row_len] = results_stm['states'][0,stm_row_index] 
        # Update FX        
        if palc_args == None:        
            FX = xconstraint-xdesired
        else:
            FX[:-1] = xconstraint-xdesired
            FX[-1] = np.matmul(np.transpose(xfree-palc_args['free_var_prev']), palc_args['delta_X*_prev']) - palc_args['delta_s'] 
        
        print('Iteration:',count+1,'|FX|=',np.linalg.norm(FX))
        
        count = count + 1 
    
    if count > Nmax:
        iterflag = True
        print('\nMaximum number of iterations exceeded. Recompute with smaller step size, different continuaton paramter, or recheck setup.\n')
   
    # Use targeted states to generate targeted Periodic orbit
    print(tf,sym_period_targ)
    results_stm = prop_cr3bp(mu, ic, tf*1/sym_period_targ, stm_bool=1, tol=int_tol, xcross_cond=0)
    
    # Compute Jacobian if retargeted orbit cannot meet while loop condition, priimarily used for Pseudo-Arc Length Continuation
    if count == 0:    
        states_final = results_stm["states"][-1,:]
        Ux, Uy, Uz, ax, ay, az = ui_partials_acc_cr3bp(mu, states_final)
        Dot_states = np.array([states_final[3], states_final[4], states_final[5], ax, ay, az])
        # Setup DF, Jacobian Matrix
        stm_temp = results_stm["stm"][stm_row_index, :, -1]
        DF[:stm_row_len,:stm_col_len] = stm_temp[:, stm_col_index]          
        if 't' in free_vars:
            DF[:stm_row_len,-1] = Dot_states[stm_row_index] 
        if sym_period_targ == 1:
            DF[:stm_row_len,:stm_col_len] = DF[:stm_row_len,:stm_col_len] - identity_temp           

    # Used for PALC
    results_stm['DF'] = DF
    results_stm['free_vars_targeted'] = xfree     
    
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

def map_vars_index_cr3bp(var_names=None):
    """
    Map variables defined in free variables or constraints to numerical index
    Dhruv Jain, Feb 25 2022
    
    Parameters
    ----------
    vars_name : list of strings
        Contains strings expressing parameters
    
    Returns
    -------
    vars_index : list of int
        Integer index code to represent the parameters
    """    
    
    variable_dict = {'x': 0, 'y': 1, 'z': 2, 'vx': 3, 'vy': 4, 'vz': 5, 'jc': 6, 't':7}
    vars_index = []
    
    if var_names != None:
        for i in range(len(var_names)):
            # Check and handle KeyError
            try: 
                vars_index.append(variable_dict[var_names[i]])
            except KeyError:
                print(vars_index[i],'is not a valid free variable')
                vars_index = 0
                break

    # Sort the indices to handle variables passed in any order
    vars_index = sorted(vars_index)
    
    return vars_index