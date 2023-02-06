"""
Created on Mon Feb 21 20:37:16 2022

@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University
    dhruvj9922@gmail.com

Objective: This file contains a Class that inherits most of the properties of the periodic_orbits class
            It is uses numerical continuation to compute a family of periodic orbit

    Features:
        1. Natural Parameter Continuation (npc) (states, time, JC)
        2. Pseduo-Arc Length Continuation (palc)

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
import pickle as pickle

import numpy as np
import scipy as sci
from cr3bp_PO_master import periodic_orbit


# Inherit attributes and methods of periodic_orbits
class periodic_orbit_fam_continuation(periodic_orbit):
    """Child class of periodic_orbot class; focuses on numerical continuation to compute family of Periodic orbits"""

    def __init__(
        self,
        sys_chars_vals,
        ic,
        tf=0,
        teval=None,
        stm_bool=0,
        xcross_cond=0,
        int_method="DOP853",
    ):
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
        # Call __init__ of cr3bp_model and periodic_orbit
        super().__init__(
            sys_chars_vals,
            ic,
            tf=tf,
            teval=teval,
            stm_bool=stm_bool,
            int_method=int_method,
        )

        # Initializes object attributes
        self.targeted_po_fam = (
            []
        )  # Stores the targeted results of a periodic orbit as a list element
        self.targeted_po_char = {
            "ic": [],
            "tf": [],
            "jc": [],
            "eigenvalues": [],
            "eigenvectors": [],
            "monodromy": [],
        }  # Stores key characterisitcs of each targeter periodic orbit in the family
        self.targeted_po_char["sys_val"] = sys_chars_vals

    def npc_po_fam(
        self,
        free_vars,
        constraints,
        sym_period_targ=1 / 2,
        conv_tol=1e-12,
        int_tol=1e-12,
        Nmax=50,
        step_size=1e-4,
        num_fam_members=1,
        param_continue="x",
        line_search=False,
        line_search_params=None,
    ):
        """
        Use Natural Parameter Continuation to compute family of periodic orbit
        => Saves data to targeted_po_char and targeted_po_fam

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
        conv_tol : TYPE, optional
            DESCRIPTION. The default is 1e-12.
        int_tol : TYPE, optional
            DESCRIPTION. The default is 1e-12.
        Nmax : TYPE, optional
            DESCRIPTION. The default is 50.
        step_size : float, optional
            Initial step size by which the parameter to be continued in will be updated to comptue a new family member. The default is 1e-4.
        num_fam_members : int, optional
            Numver of family members to be computed. The default is 1.
        param_continue : string, optional
            Parameter to continue along: usually state or time. The default is 'x'.
        line_search : Boolean, optional
            True: Update step size by a factor if unable solution not converged in set #iterations. The default is False.
        line_search_params : dict, optional
            dictionary to store parameters for line search
            'attenuation_factor' = factor by which to decrease the step size when targeter didnt converge
            'step_size0' = initial step_size
            'lower_lim_factor' = fraction of the initial step size till which the step size will be decreased
            The default is None.
        """
        self.conv_tol = conv_tol

        # Note: To target a specific JC, set periodic_orbit.obeject's JCd to be of that value

        # Perform logical checks and update free_vars & constraint based  on param_continue
        (
            check_logic_val,
            free_vars,
            constraints,
        ) = self.npc_po_logic_check_continueparam_update(
            free_vars, constraints, param_continue
        )
        if check_logic_val is None:
            return None, None

        iterflag = False
        count_fam_member = 0

        # Setup line_search_params
        if line_search_params is None:
            line_search_params = {}
            line_search_params["attenuation_factor"] = 0.8
            line_search_params["step_size0"] = step_size
            line_search_params["lower_lim_factor"] = 0.1

        # While loop to compute family members
        while count_fam_member < num_fam_members and iterflag is False:
            results, iterflag = self.single_shooter(
                free_vars,
                constraints,
                sym_period_targ=sym_period_targ,
                palc_args=None,
                conv_tol=conv_tol,
                int_tol=int_tol,
                Nmax=Nmax,
            )

            # Targeter failed to converge, use line search to update step size or end the loop
            if iterflag is True and line_search is True:
                # Use Line Search: Update Step size and reset recompute
                iterflag, step_size = self.line_search(
                    results, step_size, param_continue, line_search_params
                )

            # Targeter converged => save data and update the setup to compute the next family member
            elif iterflag is False:
                print("# PO family member = ", count_fam_member + 1, "\n")
                self.targeted_po_fam.append(
                    results
                )  # tf is updated at the end of the shooter function
                self.ic = copy.deepcopy(
                    results["states"][0, :]
                )  # To not update save data as values are passed as object reference
                self.save_targeted_po_char(results)
                count_fam_member += 1

                if count_fam_member < num_fam_members:
                    self.npc_po_update(
                        results, step_size, param_continue
                    )  # Continue the param_continue with the current step size

            else:
                print("Recheck targeter setup")
                break

    def npc_po_logic_check_continueparam_update(
        self, free_vars, constraints, param_continue
    ):
        """
        Perform logical checks for Natural Parameter Continuation function and update free_vars and constraint based on param_continue

        Parameters
        ----------
        free_vars : list of strings
            Describes the parameters set as free variables
            Possible paramters: ['x','y','z','vx','vy','vz','t']
        constraints : list of strings
            Describes the parameters set as constraints
            Possible paramters: ['x','y','z','vx','vy','vz','jc']
        param_continue : string, optional
            Parameter to continue along: usually state or time. The default is 'x'.

        Returns
        -------
        value to see if logical checks are passed:
            0-> YES,
            None -> NO
        free_vars : list of strings
            Describes the parameters set as free variables
            Possible paramters: ['x','y','z','vx','vy','vz','t']
        constraints : list of strings
            Describes the parameters set as constraints
            Possible paramters: ['x','y','z','vx','vy','vz','jc']
        """
        # Check if paramter to be continued in is defined as a free variable
        if param_continue in free_vars or param_continue == "jc":
            free_vars = copy.deepcopy(free_vars)
            constraints = copy.deepcopy(constraints)

            # Remove paramter to be continued in from free variable, if not 'jc'
            if param_continue != "jc":
                free_vars.remove(param_continue)
            # To add 'jc' to constraints if family to be continued; defined as param_continue but explicitly as a constraint
            if param_continue == "jc" and "jc" not in constraints:
                constraints.append("jc")

        else:
            print(
                "Paramter that is to be continued in is not defined as a free variable or constraint. Make sure to that the parameter to be continued in can be varied and included in free_vars/constraints"
            )
            return None, free_vars, constraints

        # Assign a value of JCd if continuing in JC but JCd is not given
        if "jc" in constraints and self.JCd is None:
            self.JCd = self.JC(self.ic)
        print("JCd", self.JCd)
        return 0, free_vars, constraints

    def npc_po_update(self, results, step_size, param_continue):
        """
        Update parameter to be continued in with the current step size
        Parameters
        ----------
        results: Dictionary
            't': time history, [nd]
            'states': state history, [:,i] => 6 states @ t = t_i, [nd]
            'yevents': states at prescribed event, [nd]
            'tevents': time stamp of yevents, [nd]
            'stm': STM element history, [:,:,i] => STM @ t = t_i
        step_size : float, optional
            Initial step size by which the parameter to be continued in will be updated to comptue a new family member. The default is 1e-4.
        param_continue : string
            Parameter to continue along: usually state or time.
        """
        # Map param_continue (str) to its set index
        param_conti_index = self.map_vars_index_cr3bp([param_continue])[0]

        if param_conti_index < 6:
            self.ic[param_conti_index] += step_size  # Update state
        elif param_continue == "t":
            self.tf += step_size  # Update time
        elif param_continue == "jc":
            self.JCd += step_size  # Update JC

    def palc_po_fam(
        self,
        free_vars,
        constraints,
        sym_period_targ=1 / 2,
        conv_tol=1e-12,
        int_tol=1e-12,
        Nmax=50,
        step_size=1e-4,
        num_fam_members=1,
        param_continue="x",
        line_search=False,
        line_search_params=None,
    ):
        """
        Use Pseudo-arc length Continuation to compute family of periodic orbit
        => Saves data to targeted_po_char and targeted_po_fam

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
        conv_tol : TYPE, optional
            DESCRIPTION. The default is 1e-12.
        int_tol : TYPE, optional
            DESCRIPTION. The default is 1e-12.
        Nmax : TYPE, optional
            DESCRIPTION. The default is 50.
        step_size : float, optional
            Initial step size by which the parameter to be continued in will be updated to comptue a new family member. The default is 1e-4.
        num_fam_members : int, optional
            Numver of family members to be computed. The default is 1.
        param_continue : string, optional
            Parameter to continue along: usually state or time. The default is 'x'.
        line_search : Boolean, optional
            True: Update step size by a factor if unable solution not converged in set #iterations. The default is False.
        line_search_params : dict, optional
            dictionary to store parameters for line search
            'attenuation_factor' = factor by which to decrease the step size when targeter didnt converge
            'step_size0' = initial step_size
            'lower_lim_factor' = fraction of the initial step size till which the step size will be decreased
            The default is None.
        """

        check_logic_val = self.palc_po_logic_check(
            free_vars, constraints, sym_period_targ
        )
        if check_logic_val is None:
            return None, None

        null_vec = np.ones(len(free_vars))

        # Setup PALC arguments to target PALC based orbits
        palc_args = {}
        palc_args["delta_s"] = step_size

        # Retarget as one of the parameters would have been removed from NPC or any other targetere setup and DF will not be large enough to be fully determined with PALC constraint
        retargeted_orbit, iterflag = self.single_shooter(
            free_vars,
            constraints,
            sym_period_targ=sym_period_targ,
            conv_tol=conv_tol,
            int_tol=int_tol,
            Nmax=Nmax,
        )

        if sym_period_targ == 1:
            # Setup Phase Condition, all states are free variables
            palc_args = self.palc_po_phase_constraint(
                free_vars, retargeted_orbit, palc_args
            )

            palc_args["prev_conv_soln"] = retargeted_orbit["states"][0, :]
            # Assuming 7 free var, 6 states + time
            DF = np.zeros((len(free_vars) - 1, len(free_vars)))
            DF[:-1, :] = retargeted_orbit["DF"]

            DF[-1, :-1] = palc_args[
                "dx/dtheta"
            ]  # Time phase constraint part is 0
            retargeted_orbit["DF"] = copy.deepcopy(DF)

        # Compute null vector and palc constraint components
        palc_args, null_vec = self.palc_null_vect_update(
            free_vars, retargeted_orbit, null_vec, palc_args
        )

        # Setup line_search_params
        if line_search_params is None:
            line_search_params = {}
            line_search_params["attenuation_factor"] = 0.8
            line_search_params["step_size0"] = step_size
            line_search_params["lower_lim_factor"] = 0.1

        iterflag = False
        count_fam_member = 0
        while count_fam_member < num_fam_members and iterflag is False:
            results, iterflag = self.single_shooter(
                free_vars,
                constraints,
                sym_period_targ=sym_period_targ,
                palc_args=palc_args,
                conv_tol=conv_tol,
                int_tol=int_tol,
                Nmax=Nmax,
            )

            if iterflag is True and line_search is True:
                # Use Line Search: Update Step size and reset recompute
                iterflag, step_size = self.line_search(
                    results,
                    step_size,
                    param_continue,
                    line_search_params,
                    palc_args=palc_args,
                )

            elif iterflag is False:
                print("# PO family member = ", count_fam_member + 1, "\n")
                self.targeted_po_fam.append(
                    results
                )  # tf is updated at the end of the shooter function
                self.ic = copy.deepcopy(
                    results["states"][0, :]
                )  # To not update save data as values are passed as object reference
                palc_args, null_vec = self.palc_null_vect_update(
                    free_vars, results, null_vec, palc_args
                )  # Compute null vector and palc constraint components
                palc_args = self.palc_po_phase_constraint(
                    free_vars, results, palc_args
                )  # Compute Phase constraint
                self.save_targeted_po_char(results)

                count_fam_member += 1

            else:
                print("Recheck targeter setup")
                break

        return self.targeted_po_fam, self.targeted_po_char

    def palc_po_logic_check(self, free_vars, constraints, sym_period_targ):
        """
        Performs logical check for PALC setup
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

        Returns
        -------
        value to see if logical checks are passed:
            0-> YES,
            None -> NO

        """

        print(
            "\nAssumes the details of the orbit passed are that of a targeted Periodic Orbit\n"
        )
        if "jc" in constraints:
            print("JC cannot be constrained when using PALC")
            return None

        if sym_period_targ == 1:
            null_vect_dim_check = 2
        else:
            null_vect_dim_check = 1
        if len(free_vars) != len(constraints) + null_vect_dim_check:
            print(
                "Recheck Free variable and constraint setup as Null space needs to be exactly one"
            )
            return None

        return 0

    def palc_null_vect_update(self, free_vars, results, null_vec, palc_args):
        """
        Computes the Null vector and other components for PALC constraint
        Parameters
        ----------
        free_vars : list of strings
            Describes the parameters set as free variables
            Possible paramters: ['x','y','z','vx','vy','vz','t']
        results: Dictionary
            't': time history, [nd]
            'states': state history, [:,i] => 6 states @ t = t_i, [nd]
            'yevents': states at prescribed event, [nd]
            'tevents': time stamp of yevents, [nd]
            'stm': STM element history, [:,:,i] => STM @ t = t_i
        null_vec : TYPE
            DESCRIPTION.
        palc_args : dictionary containting information for Pseudo-Arc Length Continuation, optional
            'delta_s': step-size
            'free_var_prev': Free variables of the previously converged orbit
            'delta_X*_prev': Null-Vector of DF of previously converged  orbit
            'prev_conv_soln': targeted IC of the previously converged orbit
            'dx/dtheta': phase constraint
            The default is None

        Returns
        -------
        palc_args : dictionary containting information for Pseudo-Arc Length Continuation, optional
            'delta_s': step-size
            'free_var_prev': Free variables of the previously converged orbit
            'delta_X*_prev': Null-Vector of DF of previously converged  orbit
            'prev_conv_soln': targeted IC of the previously converged orbit
            'dx/dtheta': phase constraint
            The default is None
        null_vec : numpy ndarray (len(xfree)x1)
            Null vector of DF => Tells us the tangent space of DF to get other family member

        """
        # Compute Null Space
        free_var_prev_null_vect = sci.linalg.null_space(results["DF"])
        if np.size(free_var_prev_null_vect, 1) != 1:
            print(
                "Null space is not one, nullity is",
                np.size(free_var_prev_null_vect, 1),
                "continuing with first null vector",
            )
            free_var_prev_null_vect = free_var_prev_null_vect[:, 0]

        free_var_prev_null_vect = free_var_prev_null_vect.flatten()

        # Check if sign of null vector is same as previous null vector, if not then change the sign
        null_vecs_dot = np.dot(free_var_prev_null_vect, null_vec)
        null_vec = free_var_prev_null_vect * np.sign(null_vecs_dot)
        palc_args["free_var_prev"] = results["free_vars_targeted"]
        palc_args["prev_conv_soln"] = results["states"][0, :]
        palc_args["delta_X*_prev"] = null_vec

        return palc_args, null_vec

    def palc_po_phase_constraint(self, free_vars, results, palc_args):
        """
        Computes the phase constraint for a periodic orbit

        Parameters
        ----------
        free_vars : list of strings
            Describes the parameters set as free variables
            Possible paramters: ['x','y','z','vx','vy','vz','t']
        results: Dictionary
            't': time history, [nd]
            'states': state history, [:,i] => 6 states @ t = t_i, [nd]
            'yevents': states at prescribed event, [nd]
            'tevents': time stamp of yevents, [nd]
            'stm': STM element history, [:,:,i] => STM @ t = t_i
        palc_args : dictionary containting information for Pseudo-Arc Length Continuation, optional
            'delta_s': step-size
            'free_var_prev': Free variables of the previously converged orbit
            'delta_X*_prev': Null-Vector of DF of previously converged  orbit
            'prev_conv_soln': targeted IC of the previously converged orbit
            'dx/dtheta': phase constraint
            The default is None

        Returns
        -------
        palc_args : dictionary containting information for Pseudo-Arc Length Continuation, optional
            'delta_s': step-size
            'free_var_prev': Free variables of the previously converged orbit
            'delta_X*_prev': Null-Vector of DF of previously converged  orbit
            'prev_conv_soln': targeted IC of the previously converged orbit
            'dx/dtheta': phase constraint
            The default is None
        """
        # Phase constraint = dxi/dtheta  = dxi/dt*(Period/2pi)
        _, _, _, ax, ay, az = self.ui_partials_acc_cr3bp(
            results["states"][0, :]
        )
        palc_args["dx/dtheta"] = (
            np.array(
                [
                    results["states"][0, 3],
                    results["states"][0, 4],
                    results["states"][0, 5],
                    ax,
                    ay,
                    az,
                ]
            )
            * results["t"][-1]
            / (2 * np.pi)
        )
        free_vars_index = self.map_vars_index_cr3bp(free_vars)
        stm_col_index = [
            free_vars_index[i]
            for i in range(len(free_vars))
            if free_vars_index[i] < 6
        ]
        palc_args["dx/dtheta"] = palc_args["dx/dtheta"][stm_col_index]

        return palc_args

    def line_search(
        self,
        results,
        step_size,
        param_continue,
        line_search_params,
        palc_args=None,
    ):
        """
        Updates the step size based on the given attentuation factor, till the step size is bigger than a given fraction of the initial step size

        Parameters
        ----------
        results: Dictionary
            't': time history, [nd]
            'states': state history, [:,i] => 6 states @ t = t_i, [nd]
            'yevents': states at prescribed event, [nd]
            'tevents': time stamp of yevents, [nd]
            'stm': STM element history, [:,:,i] => STM @ t = t_i
        step_size : float, optional
            Initial step size by which the parameter to be continued in will be updated to comptue a new family member. The default is 1e-4.
        param_continue : string
            Parameter to continue along: usually state or time
        line_search_params : dict
            dictionary to store parameters for line search
            'attenuation_factor' = factor by which to decrease the step size when targeter didnt converge
            'step_size0' = initial step_size
            'lower_lim_factor' = fraction of the initial step size till which the step size will be decreased
        palc_args : dictionary containting information for Pseudo-Arc Length Continuation, optional
            'delta_s': step-size
            'free_var_prev': Free variables of the previously converged orbit
            'delta_X*_prev': Null-Vector of DF of previously converged  orbit
            'prev_conv_soln': targeted IC of the previously converged orbit
            'dx/dtheta': phase constraint
            The default is None

        Returns
        -------
        iterflag: Boolean/None
            True: Targeter unable to converge
            False: Succesfully targeted
        step_size : float
            Updated value of step size that was passed as an argument to the function

        """

        # Reset continuation parameter to be the last converged value, NEGATIVE STEP SIZE
        if palc_args is None:  # Works only for NPC
            self.npc_po_update(results, -step_size, param_continue)

        # Update step size
        step_size = step_size * line_search_params["attenuation_factor"]
        print("Line search is used to update step size to:", step_size, "\n")

        # Update parameter with the new step size
        if palc_args is None:  # For NPC
            self.npc_po_update(results, step_size, param_continue)
        else:  # For PALC
            palc_args["delta_s"] = step_size

        # Check if step size is bigger than a given fraction of the initial step size
        if abs(step_size) < abs(
            line_search_params["step_size0"]
            * line_search_params["lower_lim_factor"]
        ):
            print(
                "Updated step size is too small compared to given step size. Rerun with smaller step size, attenuation factor or allowable lower limit"
            )
            iterflag = True  # STOP the WHILE loop
        else:
            iterflag = False

        return iterflag, step_size

    def save_targeted_po_char(self, results):
        """
        Appends the key characterisitcs of a newly targeted periodic orbit family member to "targeted_po_char"

        Parameters
        ----------
        results: Dictionary
            't': time history, [nd]
            'states': state history, [:,i] => 6 states @ t = t_i, [nd]
            'yevents': states at prescribed event, [nd]
            'tevents': time stamp of yevents, [nd]
            'stm': STM element history, [:,:,i] => STM @ t = t_i
        """

        # Save key characterisitcs
        self.targeted_po_char["conv_tol"] = self.conv_tol
        self.targeted_po_char["int_tol"] = self.int_tol

        self.targeted_po_char["ic"].append(
            copy.deepcopy(results["states"][0, :])
        )
        self.targeted_po_char["tf"].append(copy.deepcopy(results["t"][-1]))
        self.targeted_po_char["jc"].append(
            copy.deepcopy(self.JC(results["states"][0, :]))
        )
        self.targeted_po_char["monodromy"].append(results["stm"][:, :, -1])
        eigenvals, eigenvects = np.linalg.eig(results["stm"][:, :, -1])
        self.targeted_po_char["eigenvalues"].append(eigenvals)
        self.targeted_po_char["eigenvectors"].append(eigenvects)

    def pickle_orbit_values(self, filename):
        po_data = self.targeted_po_fam
        with open(filename + ".pickle", "wb") as f:
            pickle.dump(po_data, f, pickle.HIGHEST_PROTOCOL)

    def pickle_orbit_char(self, filename):
        po_char = self.targeted_po_char
        with open(filename + ".pickle", "wb") as f:
            pickle.dump(po_char, f, pickle.HIGHEST_PROTOCOL)
