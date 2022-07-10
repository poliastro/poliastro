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
from cr3bp_master import ui_partials_acc_cr3bp
from cr3bp_PO_targeter import map_vars_index_cr3bp


def npc_po_fam_cr3bp(
    mu,
    shooter_func,
    initial_guess,
    tf_guess,
    free_vars,
    constraints,
    sym_period_targ=1 / 2,
    JCd=None,
    conv_tol=1e-12,
    int_tol=1e-12,
    Nmax=10,
    step_size=1e-4,
    num_fam_members=1,
    param_continue="x",
    line_search=False,
):
    """
    Natural Parameter Continuation to compute a family of Periodic Orbits
    Dhruv Jain, Feb 28 2022

    Parameters
    ----------
    mu :  float, M2/(M1+M2)
        M1 and M2 are mass of Primary Bodies and M2<M1
    shooter_func : function
        Function used to target a periodic orbit
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
        Initial JC value to be targeted if param_conti is 'jc'. The default is None.
    conv_tol : float, optional
        Convergence tolerance of the constraint vector -> Acceptable tolerance of L2 norm of the constraint vector . The default is 1e-12.
    int_tol : float, optional
        Absolute = Relative Integration Tolerance
        The default is 1e-12.
    Nmax : int, optional
        Max allowable iteration of targeter while loop. The default is 50.
    step_size : float, optional
        Initial step size by which the parameter to be continued in will be updated to comptue a new family member. The default is 1e-4.
    num_fam_members : int, optional
        Numver of family members to be computed. The default is 1.
    param_continue : string, optional
        Parameter to continue along: usually state or time. The default is 'x'.
    line_search : Boolean, optional
        True: Update step size by a factor if unable solution not converged in set #iterations. The default is False.

    Returns
    -------
    targeted_po_fam: list of dictionary
         Contains list of periodic orbit dictionary elements
    targeted_po_char: dictionary
        Contains key characterisitcs of each member of the family computed:
        'ic': Targeted Initial Condition
        'tf': Period
        'jc': JC
        'monodromy': Monodromy Matrix
        'eigenvalues': Eigenvalues
        'eigenvectors': Eigenvectors
    """

    # Check if paramter to be continued in is defined as a free variable
    if param_continue in free_vars or param_continue == "jc":
        free_vars = copy.copy(free_vars)
        constraints = copy.copy(constraints)

        # Remove paramter to be continued in from free variable, if not 'jc'
        if param_continue != "jc":
            free_vars.remove(param_continue)

        # To add 'jc' to constraints if family to be continued in JC without explicitly defining it
        if param_continue == "jc" and "jc" not in constraints:
            constraints.append("jc")

        param_conti_index = map_vars_index_cr3bp([param_continue])[0]

    else:
        print(
            "Paramter that is to be continued in is not defined as a free variable or constraint. Make sure to that the parameter to be continued in can be varied and included in free_vars/constraints"
        )
        return None, None

    # Assign a value of JCd if continuing in JC but JCd is not given
    if "jc" in constraints and JCd is None:
        JCd = JC(mu, initial_guess[0:3], initial_guess[3:6])

    print("JCd", JCd)

    targeted_po_fam = []
    targeted_po_char = {
        "ic": [],
        "tf": [],
        "jc": [],
        "eigenvalues": [],
        "eigenvectors:": [],
        "monodromy": [],
    }

    iterflag = False
    count_fam_member = 0
    step_size0 = step_size

    while count_fam_member < num_fam_members and iterflag is False:

        results, iterflag = shooter_func(
            mu,
            initial_guess,
            tf_guess,
            free_vars,
            constraints,
            sym_period_targ=sym_period_targ,
            JCd=JCd,
            conv_tol=conv_tol,
            int_tol=int_tol,
            Nmax=Nmax,
        )

        if iterflag is True:
            # Use Line Search: Update Step size and recompute
            if line_search is True:
                step_size = step_size * 0.8
                print(
                    "Line search is used to update step size to:",
                    step_size,
                    "\n",
                )
                if param_conti_index < 6:
                    initial_guess[param_conti_index] -= step_size
                elif param_continue == "t":
                    tf_guess -= step_size
                elif param_continue == "jc":
                    JCd -= step_size

                if abs(step_size) < abs(step_size0 * 0.1):
                    print(
                        "Updated step size is too small compared to given step size. Rerun with smaller step size"
                    )
                else:
                    iterflag = False
        elif iterflag is False:
            print("# PO family member = ", count_fam_member + 1, "\n")
            targeted_po_fam.append(results)
            tf_guess = results["t"][-1]
            initial_guess = copy.copy(
                results["states"][0, :]
            )  # To not update save data as values are passed as object reference
            if param_conti_index < 6:
                initial_guess[param_conti_index] += step_size
            elif param_continue == "t":
                tf_guess += step_size
            elif param_continue == "jc":
                JCd += step_size

            # Save key characterisitcs
            targeted_po_char["ic"].append(copy.copy(results["states"][0, :]))
            targeted_po_char["tf"].append(copy.copy(results["t"][-1]))
            targeted_po_char["jc"].append(
                copy.copy(
                    JC(
                        mu,
                        results["states"][0, 0:3],
                        results["states"][0, 3:6],
                    )
                )
            )
            targeted_po_char["monodromy"].append(results["stm"][:, :, -1])
            eigenvals, eigenvects = np.linalg.eig(results["stm"][:, :, -1])
            targeted_po_char["eigenvalues"].append(eigenvals)
            targeted_po_char["eigenvectors:"].append(eigenvects)

            count_fam_member += 1
        else:
            print("Recheck targeter setup")
            break

    return targeted_po_fam, targeted_po_char


def palc_po_fam_cr3bp(
    mu,
    shooter_func,
    targeted_orbit,
    free_vars,
    constraints,
    sym_period_targ=1 / 2,
    conv_tol=1e-12,
    int_tol=1e-12,
    Nmax=10,
    step_size=1e-4,
    num_fam_members=1,
    line_search=False,
):
    """
    Pseudo-Arc Length Continuation to compute a family of Periodic Orbits
    Dhruv Jain, Feb 28 2022

    Parameters
    ----------
    mu :  float, M2/(M1+M2)
        M1 and M2 are mass of Primary Bodies and M2<M1
    shooter_func : function
        Function used to target a periodic orbit
    targeted_orbit : dictionary
        contains states, time and stm information to retarget the orbit to get [DF]
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
    conv_tol : float, optional
        Convergence tolerance of the constraint vector -> Acceptable tolerance of L2 norm of the constraint vector . The default is 1e-12.
    int_tol : float, optional
        Absolute = Relative Integration Tolerance
        The default is 1e-12.
    Nmax : int, optional
        Max allowable iteration of targeter while loop. The default is 50.
    step_size : float, optional
        Initial step size by which the parameter to be continued in will be updated to comptue a new family member. The default is 1e-4.
    num_fam_members : int, optional
        Numver of family members to be computed. The default is 1.
    line_search : Boolean, optional
        True: Update step size by a factor if unable solution not converged in set #iterations. The default is False.

    Returns
    -------
    targeted_po_fam: list of dictionary
         Contains list of periodic orbit dictionary elements
    targeted_po_char: dictionary
        Contains key characterisitcs of each member of the family computed:
        'ic': Targeted Initial Condition
        'tf': Period
        'jc': JC
        'monodromy': Monodromy Matrix
        'eigenvalues': Eigenvalues
        'eigenvectors': Eigenvectors

    """
    print(
        "\nAssumes the details of the orbit passed are that of a targeted Periodic Orbit\n"
    )

    if "jc" in constraints:
        print("JC cannot be constrained when using PALC")
        return None, None

    if sym_period_targ == 1:
        null_vect_dim_check = 2
    else:
        null_vect_dim_check = 1
    if len(free_vars) != len(constraints) + null_vect_dim_check:
        print(
            "Recheck Free variable and constraint setup as Null space needs to be exactly one"
        )
        return None, None

    # Retarget as one of the parameters would have been removed from NPC or any other targetere setup and DF will not be large enough to be fully determined with PALC constraint
    retargeted_orbit, iterflag = shooter_func(
        mu,
        targeted_orbit["states"][0, :],
        targeted_orbit["t"][-1],
        free_vars,
        constraints,
        sym_period_targ=sym_period_targ,
    )
    null_vec = np.ones(len(free_vars))

    # Setup PALC arguments to target PALC based orbits
    palc_args = {}
    palc_args["delta_s"] = step_size

    if sym_period_targ == 1:
        # Setup Phase Condition, all states are free variables
        _, _, _, ax, ay, az = ui_partials_acc_cr3bp(
            mu, retargeted_orbit["states"][0, :]
        )
        palc_args["dx/dtheta"] = (
            np.array(
                [
                    retargeted_orbit["states"][0, 3],
                    retargeted_orbit["states"][0, 4],
                    retargeted_orbit["states"][0, 5],
                    ax,
                    ay,
                    az,
                ]
            )
            * retargeted_orbit["t"][-1]
            / (2 * np.pi)
        )
        free_vars_index = map_vars_index_cr3bp(free_vars)
        stm_col_index = [
            free_vars_index[i]
            for i in range(len(free_vars))
            if free_vars_index[i] < 6
        ]
        palc_args["dx/dtheta"] = palc_args["dx/dtheta"][stm_col_index]

        palc_args["prev_conv_soln"] = retargeted_orbit["states"][0, :]
        # Assuming 7 free var, 6 states + time
        DF = np.zeros((len(free_vars) - 1, len(free_vars)))
        DF[:-1, :] = retargeted_orbit["DF"]

        DF[-1, :-1] = palc_args["dx/dtheta"]  # Time phase constraint part is 0
        retargeted_orbit["DF"] = copy.copy(DF)

    # Compute Null Space
    free_var_prev_null_vect = sci.linalg.null_space(retargeted_orbit["DF"])
    if np.size(free_var_prev_null_vect, 1) != 1:
        print(
            "Null space is not one, nullity is",
            np.size(free_var_prev_null_vect, 1),
            "continuing with first null vector",
        )
        free_var_prev_null_vect = free_var_prev_null_vect[0]

    free_var_prev_null_vect = free_var_prev_null_vect.flatten()

    # Check if sign of null vector is same as previous null vector, if not then change the sign
    null_vecs_dot = np.dot(free_var_prev_null_vect, null_vec)
    null_vec = free_var_prev_null_vect * np.sign(null_vecs_dot)
    palc_args["free_var_prev"] = retargeted_orbit["free_vars_targeted"]
    palc_args["delta_X*_prev"] = null_vec

    targeted_po_fam = []
    targeted_po_char = {
        "ic": [],
        "tf": [],
        "jc": [],
        "eigenvalues": [],
        "eigenvectors:": [],
        "monodromy": [],
    }

    iterflag = False
    count_fam_member = 0
    step_size0 = step_size
    initial_guess = retargeted_orbit["states"][0, :]
    tf_guess = retargeted_orbit["t"][-1]

    while count_fam_member < num_fam_members and iterflag is False:

        results, iterflag = shooter_func(
            mu,
            initial_guess,
            tf_guess,
            free_vars,
            constraints,
            sym_period_targ=sym_period_targ,
            palc_args=palc_args,
            conv_tol=conv_tol,
            int_tol=int_tol,
            Nmax=Nmax,
        )

        if iterflag is True:
            # Use Line Search: Update Step size and recompute
            if line_search is True:
                step_size = step_size * 0.8
                print(
                    "Line search is used to update step size to:",
                    step_size,
                    "\n",
                )
                palc_args["delta_s"] = step_size

                if abs(step_size) < abs(step_size0 * 0.1):
                    print(
                        "Updated step size is too small compared to given step size. Rerun with smaller step size"
                    )
                else:
                    iterflag = False
        elif iterflag is False:
            print("# PO family member = ", count_fam_member + 1, "\n")
            targeted_po_fam.append(results)
            tf_guess = results["t"][-1]
            initial_guess = copy.copy(
                results["states"][0, :]
            )  # To not update save data as values are passed as object reference

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

            _, _, _, ax, ay, az = ui_partials_acc_cr3bp(
                mu, results["states"][0, :]
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
            free_vars_index = map_vars_index_cr3bp(free_vars)
            stm_col_index = [
                free_vars_index[i]
                for i in range(len(free_vars))
                if free_vars_index[i] < 6
            ]
            palc_args["dx/dtheta"] = palc_args["dx/dtheta"][stm_col_index]

            # Save key characterisitcs
            targeted_po_char["ic"].append(copy.copy(results["states"][0, :]))
            targeted_po_char["tf"].append(copy.copy(results["t"][-1]))
            targeted_po_char["jc"].append(
                copy.copy(
                    JC(
                        mu,
                        results["states"][0, 0:3],
                        results["states"][0, 3:6],
                    )
                )
            )
            targeted_po_char["monodromy"].append(results["stm"][:, :, -1])
            eigenvals, eigenvects = np.linalg.eig(results["stm"][:, :, -1])
            targeted_po_char["eigenvalues"].append(eigenvals)
            targeted_po_char["eigenvectors:"].append(eigenvects)

            count_fam_member += 1

        else:
            print("Recheck targeter setup")
            break

    return targeted_po_fam, targeted_po_char
