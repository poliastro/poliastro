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
import copy

import numpy as np
from cr3bp_master import prop_cr3bp


def multi_shooter_nodes_setup_cr3bp(
    mu, ic, tf, n_node, node_place_opt="time", int_tol=1e-12
):
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
