"""
Created on 20 Mar 2022

@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University
        dhruvj9922@gmail.com

Objective: This file contains a Class that inherits most of the properties of cr3bp_model class
            It is used to compute a periodic orbit

    Features:
        1. Single Shooter Targeter:
            -> Variable time targeter capable of targeting states, JC and periodicity
            -> Pseudo-Arc Length Compatible with Phase coonstraint
        2. Newton-Raphson Solver
        3. Setup multiple nodes to use a Multiple Shooter Targeter
        4. Compute local and global manifolds
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
from cr3bp_model_master import cr3bp_model

# Inherit attributes and methods of cr3bp_model
class periodic_orbit(cr3bp_model):    
    """ Child class of cr3bp_model; focuss on periodic orbit computation
    """   
    
    def __init__(self, sys_chars_vals, ic, tf=0, teval=None, stm_bool=0, xcross_cond=0, int_method='DOP853'):
        """
        Constructor
        Parameters
        ----------
        sys_chars_vals : object
            object of Class sys_char
        ic : numpy ndarray (6x1), {Can handle all 42 states for CR3BP+STM integration}
            States are defined about the barycenter of the two primaries, P1 and P2
            Initial condition: 6 states to compute a trajectory;
            [0:x0, 1:y0, 2:z0, 3:vx0, 4:vy0, 5:vz0] [non-dimensional] [nd]
            default = [0,0,0,0,0,0]
        tf : float
            Integration time [nd]
            Can be negative or positive, negative => Integration in backwards time, default = 0
        teval: list of float
            Contains time stamps at which the numerically integrated results will be saved
        int_tol : float
            Absolute = Relative Integration Tolerance
            The default is 1e-12.
        stm_bool : boolean, optional
            0: CR3BP EOM Integration
            1: CR3BP + STM EOM Integration
            The default is 0.
        xcross_cond : int, optional
            0 => No y-crossing check
            1 => Events function to check when crossed y axis in any direction
            2 => Events function to check when corssed y axis from -y to +y
            The default is 0.
        int_method : string, optional
            Specify integration scheme: 'DOP853' or 'LSODA'
            The default is 'DOP853'.
        """        
        # Call __init__ of cr3bp_model
        super().__init__(sys_chars_vals, ic, tf=tf, teval=teval, stm_bool=stm_bool, int_method=int_method)
        
        # Initializes object attributes
        self.sys_period_targ = 1
        self.JCd = None
        self.xfree = np.zeros((0))
        self.xconstraints = np.zeros((0))
        self.xdesired = np.zeros((0))
        self.stm_row_index = 0
        self.stm_col_index = 0
        self.stm_row_len = 0
        self.stm_col_len = 0


    # Set initial guess, tf_guess as ic and tf
    def single_shooter(self,
        free_vars,
        constraints,
        sym_period_targ=1 / 2,
        palc_args=None,
        int_tol=1e-12,
        conv_tol=1e-12,
        Nmax=50,
    ):
        """Single shooter targeter for a Periodic Orbit defined in the CR3BP model
        Dhruv Jain, 26 Feb 2022
    
        Parameters
        ----------
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
            'delta_X*_prev': Null-Vector of DF of previously converged  orbit
            'prev_conv_soln': targeted IC of the previously converged orbit
            'dx/dtheta': phase constraint 
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
        results_stm: Targeted Periodic Orbit - Dictionary
            't': time history, [nd]
            'states': state history, [:,i] => 6 states @ t = t_i, [nd]
            'yevents': states at prescribed event, [nd]
            'tevents': time stamp of yevents, [nd]
            'stm': STM element history, [:,:,i] => STM @ t = t_i
            'DF': Jacobian Matrix, without PALC
            'free_vars_targeted': Targeted Free Varaible Vector
        iterflag: Boolean/None
            True: Targeter unable to converge
            None: Targeter Setup is incorrect
            False: Succesfully targeted
        """
        
        # Initialization
        self.sym_period_targ = sym_period_targ
        iterflag = None
        ic_original = copy.deepcopy(self.ic)
        tf_original = copy.deepcopy(self.tf)
        
        # Free vairables and Constraints paramters logical error check
        if "jc" in free_vars:
            print("Jacobi Constant cannot be a free variable")
            return 0, iterflag
        if "t" in constraints:
            print("Time should not be a constraint, instead make the tf_guess to be desired time and not add 't' as a free variable")
            return 0, iterflag
        if sym_period_targ not in [1 / 4, 1 / 2, 1]:
            print("Not a valid fraction of period to target")
            return 0
        if palc_args is not None:
            print("PALC constraint included")
        
        int_method_orig = self.int_method
        self.int_method = 'DOP853'
        # Computes tf using events function to check y axis crossing; used to get accurate TOF if symmetry is to be used 
        self.__events_func_ycrossing_check(free_vars)
        self.int_method = int_method_orig

        # Propagate I.C. to new tf 
        self.xcross_cond = 0
        self.stm_bool = 1    
        results_stm = self.propagate()
    
        # Setup Targeter components: X, FX, DF, DG
        self.__free_vars_setup_update(free_vars, purpose='setup') # X
        identity_temp = self.__constraints_setup_update(constraints, results_stm,purpose='setup') # xconstraint and xdesired
        FX = self.__FX_setup_update(palc_args) # FX = xconstraint - xdesired
        self.DF, DG = self.__DF_DG_setup(self.xfree, self.xconstraint, palc_args) # DF, DG
            
        print("FX:",constraints,"X:",free_vars,"\nTf0:",self.tf, "JC0:", self.JC(self.ic))    
        print("Iteration:", 0, "|FX|=", np.linalg.norm(FX))
    
        count = 0
        iterflag = False
        
        # Targeter loop: Use L2-norm of constraint vector and #iterations as stopping condition
        while np.linalg.norm(FX) > conv_tol and count <= Nmax:            
            # Update STM after each iteration to converge faster

            # Compute DF, Jacobian Matrix
            self.__compute_DF(free_vars, constraints, results_stm, identity_temp, palc_args)            
    
            # Update Free variable vector, include PALC constraint if PALC is being used
            if palc_args is None:
                self.xfree = self.newton_raphson_update(self.xfree, FX, self.DF)
            else:                       
                DG[:-1, :] = self.DF
                DG[-1, :] = palc_args["delta_X*_prev"] # Augmented Jacobian matrix with PALC constraint
                self.xfree = self.newton_raphson_update(self.xfree, FX, DG)
    
            # Update IC and tf using the updated xfree
            self.__free_vars_setup_update(free_vars)
            
            # Propagate with updated IC and tf
            results_stm = self.propagate()
            
            # Update xconstraint
            _ = self.__constraints_setup_update(constraints, results_stm)

            FX = self.__FX_setup_update(palc_args)
    
            print("Iteration:", count + 1, "|FX|=", np.linalg.norm(FX))
    
            count = count + 1
            
        # Use targeted states to generate targeted Periodic orbit
        print(self.tf, sym_period_targ)
        
        self.tf = self.tf * 1 / sym_period_targ
        results_stm = self.propagate()
        
        # Compute Jacobian if retargeted orbit cannot meet while loop condition, priimarily used for Pseudo-Arc Length Continuation
        if count == 0:
            self.__compute_DF(free_vars, constraints, results_stm, identity_temp, palc_args)       
            # dxf_dt = self.compute_dxf_dt(results_stm)    
            # # Setup DF, Jacobian Matrix
            # stm_temp = results_stm["stm"][self.stm_row_index, :, -1]
            # self.DF[:self.stm_row_len, :self.stm_col_len] = stm_temp[:, self.stm_col_index]
            # if "t" in free_vars:
            #     self.DF[:self.stm_row_len, -1] = dxf_dt[self.stm_row_index]
            # if sym_period_targ == 1:
            #     self.DF[:self.stm_row_len, :self.stm_col_len] = (self.DF[:self.stm_row_len, :self.stm_col_len] - identity_temp)
    
        if count > Nmax:
            iterflag = True
            print("\nMaximum number of iterations exceeded. Recompute with smaller step size, different continuaton paramter, or recheck setup.\n")
            self.tf = tf_original
            self.ic = ic_original
    
        # Used for PALC
        results_stm["DF"] = self.DF
        results_stm["free_vars_targeted"] = self.xfree
        
        return results_stm, iterflag
    
    def __events_func_ycrossing_check(self, free_vars):
        """
        If symmetry is to be used to target then use events function to get "tf" till next y-axis crossinng        
        Parameters
        ----------
        free_vars : list of strings
            Describes the parameters set as free variables
            Possible paramters: ['x','y','z','vx','vy','vz','t']
        """
   
        # Use events function to compute y-crossing if symmetry is to be used for targeter
        if "t" in free_vars and self.sym_period_targ != 1:
            self.xcross_cond = 1
            results_stm = self.propagate()
            if len(results_stm["tevents"][0][:]) > 1:
              # condition so that it works for corrected solutions
                self.tf = results_stm["tevents"][0][1]
            else:  # Assumes that IC doesnt count as first pass, the next crossing that is needed for targeting is the yevent
                self.tf = results_stm["tevents"][0][0]
        else:  # If need to used time fixed targeter or Periodicity targeter
            self.tf = self.tf * self.sym_period_targ
    
    
    def __free_vars_setup_update(self, free_vars, purpose='update'):
        """
        Computes intial Free Variable Vector and Updates IC and tf after computing a new xfree
        Parameters
        ----------
        free_vars : list of strings
            Describes the parameters set as free variables
            Possible paramters: ['x','y','z','vx','vy','vz','t']
        purpose : string, optional
            to check the purpose for which the function is to be used - to 'setup' or 'update'. The default is 'update'.
        """
        if purpose != 'setup' and purpose != 'update':
            print('Incorrect purpose passed')
            self.ic = 0
            self.tf = 0
        
        else:
            # Setup Free Varaiable Vector
            if purpose == 'setup':
                # Map paramter strings to index
                free_vars_index = self.map_vars_index_cr3bp(free_vars)
                self.xfree = np.zeros(len(free_vars_index))
                
                self.stm_col_index = [free_vars_index[i] for i in range(len(self.xfree)) if free_vars_index[i] < 6]
                self.stm_col_len = len(self.stm_col_index)
        
                if "t" in free_vars:
                    self.xfree[:-1] = self.ic[self.stm_col_index]
                    self.xfree[-1] = self.tf
                else:
                    self.xfree = self.ic[self.stm_col_index]
    
            # Update Free Varaiable Vector        
            elif purpose == 'update':            
                if "t" in free_vars:
                    self.ic[self.stm_col_index] = self.xfree[:-1]
                    self.tf = self.xfree[-1]
                else:
                    self.ic[self.stm_col_index] = self.xfree
            

    
    def __constraints_setup_update(self, constraints, results_stm, purpose ='update'):
        """
        Computes intial xconstraint and xdesired Vector and Updates the two after computing a new xfree
        Parameters
        ----------
        constraints : list of strings
            Describes the parameters set as constraints
            Possible paramters: ['x','y','z','vx','vy','vz','jc']
        result_stm: Dictionary
            't': time history, [nd]
            'states': state history, [:,i] => 6 states @ t = t_i, [nd]
            'yevents': states at prescribed event, [nd]
            'tevents': time stamp of yevents, [nd]
            'stm': STM element history, [:,:,i] => STM @ t = t_i
        purpose : string, optional
            to check the purpose for which the function is to be used - to 'setup' or 'update'. The default is 'update'.
        """    
        
        if purpose != 'setup' and purpose != 'update':
            print('Incorrect purpose passed')
            self.ic = 0
            self.tf = 0
            identity_temp = 0

        else:
            # Initial Setup
            if purpose == 'setup':
                # Map paramter strings to index
                constraints_index = self.map_vars_index_cr3bp(constraints)
                self.xconstraint = np.zeros(len(constraints_index))
                self.xdesired = np.zeros(len(constraints_index))
                self.stm_row_index = [constraints_index[i] for i in range(len(self.xconstraint)) if constraints_index[i] < 6]
                self.stm_row_len = len(self.stm_row_index)
    
            # Setup/Update Constraint Vector and Desired vector, xconstraint where FX = xconstraint - xdesired
            if "jc" in constraints:
                self.xconstraint[:-1] = results_stm["states"][-1, self.stm_row_index]
                self.xconstraint[-1] = self.JC(self.ic) # JC of IC
                self.xdesired[-1] = self.JCd  # Desired JC
            else:
                self.xconstraint = results_stm["states"][-1, self.stm_row_index]
        
            # Updated desired to be inital state of orbit
            if self.sym_period_targ == 1:
                self.xdesired[:self.stm_row_len] = results_stm["states"][0, self.stm_row_index]
                # Create identity like matrix to be subtracted from STM
                identity_mat = np.eye(6)
                temp = identity_mat[self.stm_row_index, :]
                identity_temp = temp[:, self.stm_col_index]
            else:
                identity_temp = 0
        
        return identity_temp
    
    
    def __FX_setup_update(self, palc_args):
        """
        Computes the Constraint vector = F(x) = xconstraint-xdesired
        Parameters
        ----------
        palc_args : dictionary containting information for Pseudo-Arc Length Continuation, optional
            'delta_s': step-size
            'free_var_prev': Free variables of the previously converged orbit
            'delta_X*_prev': Null-Vector of DF of previously converged  orbit
            'prev_conv_soln': targeted IC of the previously converged orbit
            'dx/dtheta': phase constraint 
            The default is None

        Returns
        -------
        FX : numpy ndarray (length_constraints + 1 or 2 x 1)
            Constraint Vector, constraints of all states, includes phase and PALC constraint
        """ 
        # Create FX
        if palc_args is None:
            FX = self.xconstraint - self.xdesired
        else:
            # Check if need to add phase constraint
            if self.sym_period_targ == 1:
                FX = np.zeros((len(self.xconstraint) + 2))
                FX[:-2] = self.xconstraint - self.xdesired
                
                # Phase Constraint
                FX[-2] = np.dot(self.ic[self.stm_col_index]-palc_args['prev_conv_soln'][self.stm_col_index],palc_args['dx/dtheta'])
                
            else:
                FX = np.zeros((len(self.xconstraint) + 1))
                FX[:-1] = self.xconstraint - self.xdesired
                
            # PALC Constraint
            FX[-1] = (np.matmul(np.transpose(self.xfree - palc_args["free_var_prev"]), palc_args["delta_X*_prev"]) - palc_args["delta_s"])

        return FX


    def __DF_DG_setup(self, xfree, xconstraint, palc_args):
        if palc_args is not None and self.sym_period_targ == 1:
            DF = np.zeros((len(xconstraint)+1, len(xfree))) # To account for phase constraint when PALC used
        else:    
            DF = np.zeros((len(xconstraint), len(xfree)))
        DG = np.zeros((len(xfree), len(xfree)))  # DF for PALC
    
        return DF, DG

    def compute_dxf_dt(self, results_stm): 
        """
        Computes the time derivative of states @ t = tf 
        Parameters
        ----------
        result_stm: Dictionary
            't': time history, [nd]
            'states': state history, [:,i] => 6 states @ t = t_i, [nd]
            'yevents': states at prescribed event, [nd]
            'tevents': time stamp of yevents, [nd]
            'stm': STM element history, [:,:,i] => STM @ t = t_i
        Returns
        -------
        dxf_dt : numpy ndarray (6x1)
            d[xf, yf, zf, vxf, vyf, vzf]/dt
        """
        
        states_final = results_stm["states"][-1, :]
        Ux, Uy, Uz, ax, ay, az = self.ui_partials_acc_cr3bp(states_final)
        dxf_dt = np.array([states_final[3], states_final[4], states_final[5], ax, ay, az])
        return dxf_dt

    def compute_dJC_dxi(self):
        """
        Computes the analytical elements of Jacobian for "dJC/dstate"
        Returns
        -------
        dJC_dx : numpy ndarray (6x1), float
            dJC/d[x,y,z,vx,vy,vz]
        """
        # Compute d(JC)/d(x)
        Ux_ic, Uy_ic, Uz_ic, _, _, _ = self.ui_partials_acc_cr3bp(self.ic)
        dJC_dx = np.array([2 * Ux_ic, 2 * Uy_ic, 2 * Uz_ic, -2 * self.ic[3], -2 * self.ic[4], -2 * self.ic[5]])
        return dJC_dx

    def __compute_DF(self,free_vars,constraints,results_stm, identity_temp,palc_args):    
        """
        Computes the Analytical Jacobian Matrix, "DF", for a PO targeter
        Parameters
        ----------
        free_vars : list of strings
            Describes the parameters set as free variables
            Possible paramters: ['x','y','z','vx','vy','vz','t']
        constraints : list of strings
            Describes the parameters set as constraints
            Possible paramters: ['x','y','z','vx','vy','vz','jc']
        results_stm : TYPE
            DESCRIPTION.
        identity_temp : TYPE
            DESCRIPTION.
        palc_args : dictionary containting information for Pseudo-Arc Length Continuation, optional
            'delta_s': step-size
            'free_var_prev': Free variables of the previously converged orbit
            'delta_X*_prev': Null-Vector of DF of previously converged solution
            The default is None
        """
        dxf_dt = self.compute_dxf_dt(results_stm)    
        dJC_dx = self.compute_dJC_dxi()
        # Extract required elements of the Monodromy matrix for the targeter setup
        stm_temp = results_stm["stm"][self.stm_row_index, :, -1]

        # Setup DF, Jacobian Matrix
        self.DF[:self.stm_row_len, :self.stm_col_len] = stm_temp[:, self.stm_col_index]
        if "jc" in constraints:
            self.DF[-1, :self.stm_col_len] = dJC_dx[self.stm_col_index]
        if "t" in free_vars:
            self.DF[:self.stm_row_len, -1] = dxf_dt[self.stm_row_index]
        if self.sym_period_targ == 1:
            self.DF[:self.stm_row_len, :self.stm_col_len] = (self.DF[:self.stm_row_len, :self.stm_col_len] - identity_temp)
            # Account for Phase constraint with PALC
            if palc_args is not None:
                self.DF[-1,:len(palc_args['dx/dtheta'])] = palc_args['dx/dtheta']


    def newton_raphson_update(self, xfree, FX, DG):
        """
        Multi-dimensional Newton-Raphson Method to update inital guess
        
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
        else:  # Not a square matrix, use Pseudo-inverse
            DG_inv = np.linalg.pinv(DG)
    
        delta_x_guess = np.matmul(DG_inv, FX)
        xfree = xfree - delta_x_guess
    
        return xfree
    
    
    def map_vars_index_cr3bp(self, var_names=None):
        """
        Map variables defined in free variables or constraints to numerical index
    
        Parameters
        ----------
        vars_name : list of strings
            Contains strings expressing parameters
    
        Returns
        -------
        vars_index : list of int
            Integer index code to represent the parameters
        """
    
        variable_dict = {"x": 0, "y": 1, "z": 2, "vx": 3, "vy": 4, "vz": 5, "jc": 6, "t": 7}
        vars_index = []
    
        if var_names is not None:
            for i in range(len(var_names)):
                # Check and handle KeyError
                try:
                    vars_index.append(variable_dict[var_names[i]])
                except KeyError:
                    print(vars_index[i], "is not a valid free variable")
                    vars_index = 0
                    break
    
        # Sort the indices to handle variables passed in any order
        vars_index = sorted(vars_index)
    
        return vars_index
            
    def multi_shooter_nodes_setup_cr3bp(self, n_node, node_place_opt="time"):
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
    
    
        Parameters
        ----------
        n_node : int, Number of nodes
            IC is node 1 and nth node is node when propagated by tn value should reach the vicinity of the desired final state
        node_place_opt : string, optional
            'time': Computes the n_nodes placed after NEARLY equal time intervals
            'index': Computes the n_nodes placed after NEARLY equal time history INDEX
    
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
    
        self.stm_bool = 0
        self.xcross_cond = 0
        results_stm = self.propagate() 
        ic_node = []
        t_node = []
    
        ic_node.append(self.ic)
    
        if node_place_opt == "time":
            # Time of each segment
            ti = np.linspace(0, self.tf, n_node + 1)
            ti = ti[
                1:
            ]  # As ti is the time from node_i to next node, omit the first time, i.e. 0
    
            for node_counter in range(1, n_node):
                index = np.argmin(
                    abs(results_stm["t"] - ti[node_counter - 1])
                )  # compute index when index from time history when time is nearly
                ic_node.append(results_stm["states"][index, :])
                t_node.append(
                    results_stm["t"][index] - sum(t_node)
                )  # t_node runs one node behind as it is the time to next node, subtract sum as ti+t2+t3 = tf
    
            t_node.append(
                self.tf - sum(t_node)
            )  # Time from final node to vicinity of desired state
    
        elif node_place_opt == "index":
            num_index = len(results_stm["t"])
            indices = np.linsapce(
                0, num_index, n_node, dtype="int"
            )  # Linearly spaced indices
    
            for node_counter in range(1, n_node):
                ic_node.append(results_stm["states"][indices[node_counter], :])
                t_node.append(
                    results_stm["t"][indices[node_counter]] - sum(t_node)
                )  # t_node runs one node behind as it is the time to next node, subtract sum as ti+t2+t3 = tf
    
            t_node.append(
                self.tf - sum(t_node)
            )  # Time from final node to vicinity of desired state
    
        return ic_node, t_node
    
    def propagate_po_n_rev(self,po_ic = None, n_rev=1):
        """
        Since the results may be changed when  propagate() for t != n*tf is called, so use this func to get data after n revs
        Assumes that ic are that of targeted PO

        Parameters
        ----------
        n_rev : int, optional
            Number of revs of PO. The default is 1.
            
        Returns
        -------
        Updates self.results
        """
        if po_ic is None:
            po_ic = self.ic
      
        return self.propagate(ic=po_ic, tf=n_rev*self.tf)
    
    
    def local_manifold_gen(self, dist, eigenvect_manifold, prop_time, teval_len = 1000, po_ic = None,get_initial_states=False):
        """
        Computes local manifold of a PO
        
        Parameters
        ----------
        dist : float, km
            DESCRIPTION.
        eigenvect_manifold : numpy ndarray
            Eigenvector use to compute the manifold; size: 6 x 1
        prop_time : float
            Propagation time of the manifold. +ve for unstable, -ve for stable
        teval_len : T, optional
            Number of time steps, used to get manifold states at set time. The default is 100.
        po_ic : ndarray, optional
            PO states. Size: 6. The default is None.
        
        Returns
        -------
        local_manifold_states : numpy ndarray
            It saves the states of manfiolds of a single invariant curve
            [teval length, 6 states]
        """
        if po_ic is None:
            po_ic = self.ic
        
        dist = dist/self.sys_chars_vals.lstar# Non-dimensionalize dist from km to [nd]
        
        manifold_step_off = dist * eigenvect_manifold/np.linalg.norm(eigenvect_manifold[:3])
        x_ic_manifold = po_ic + manifold_step_off
        
        if get_initial_states is False:        
            local_manifold_states = np.zeros((teval_len,6))
            self.teval = np.linspace(0, abs(prop_time), teval_len)*np.sign(prop_time)
            
            int_method_orig = self.int_method
            self.int_method = 'DOP853'
            # self.teval is implicitly called as integrator is set to DOP853
            self.propagate(ic=x_ic_manifold, stm_bool=0,tf=prop_time)
            self.int_method = int_method_orig
        
            local_manifold_states =  copy.deepcopy(self.results['states'])
        else:
            local_manifold_states = copy.deepcopy(x_ic_manifold)
            
        return local_manifold_states
    
    def global_manifold_gen(self, dist, eigenvect_manifold, prop_time, teval_len = 100, po_ic=None, num_po_states=20, get_initial_states=False, start_time = 0, stop_time = None):
        """        
        Computes the global manifold in a certain direction (+ or -) of a certain type(stable or unstable) of a QPO
            - Change the sign of dist to get the other hald of the manifold
            - Change the eigenvect_manifold to get the other type of manifold 
                and make sure that the sign of prop_time is appropriate
        
        Parameters
        ----------
        dist, eigenvect_manifold, prop_time, teval_len, po_ic:, get_initial_states Inputs for local_manifold_gen()
        num_po_states : int, optional
            Number of PO states. It is used to compute num_po_states which are equally spaced by Time around a PO.
            The default is 20.
        start_time: float, optional
            Time at which the first manifold will be computed, useful when manifolds between set time interval is to be computed
        stop_time: float, optional
            Time at which the last manifold will be computed, useful when manifolds between set time interval is to be computed        

        Returns
        -------
        global_manifold_states: numpy ndarray
            It saves the states of manfiolds of invariant curves
            [number of invariant curve, 6 states]
        """
        if po_ic is None:
            po_ic = self.ic
        
        if stop_time is None or stop_time == 0:
            stop_time = self.tf

        tf_po_state = np.linspace(start_time, stop_time, num_po_states, endpoint = False) # Enpoint defaults to true. last value will be replaced by first to help with surface plot
        
        global_manifold_states = np.zeros((num_po_states,teval_len,6))    
        
        for count_po_states in range(num_po_states):
            
            # Obtain states at dififerent times along PO and STM(tf,0)
            int_method_orig = self.int_method
            self.int_method = 'boost'
            result_stm = self.propagate(ic=po_ic, tf=tf_po_state[count_po_states])
                    
            # eigenvect of nth invariant curve = STM(tf_i,0)*eigenvector of the first invar curve
            eigvector_tfi = np.real(np.matmul(copy.deepcopy(result_stm['stm'][:,:,-1]), eigenvect_manifold))
            
            self.int_method = int_method_orig
            
            global_manifold_states[count_po_states,:,:] = self.local_manifold_gen(dist = dist, eigenvect_manifold = eigvector_tfi, 
                                                                                  prop_time = prop_time, teval_len=teval_len, 
                                                                                  po_ic = result_stm['states'][-1,:],
                                                                                  get_initial_states=get_initial_states)
            
        return global_manifold_states
