"""
Created on Mon Feb 21 20:37:16 2022

@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University
    dhruvj9922@gmail.com

Objective: This file contains functions required to target a Periodic Orbit in the Circular Restricted
          Three Body Problem (CR3BP) model

    Features:
       1. Compute inital guess for L1, L2, L3 Lyapunov orbits using linear approximation

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
from cr3bp_lib_JC_calc import lib_pt_loc
from cr3bp_master import rel_dist_cr3bp, uii_partials_cr3bp


def ig_lyap_orb_collinear_li_cr3bp(mu, pert_x, lib_num=1):
    """
    Calculate the linear inital guess of a lyapunov orbit by calculating the
    linear approximation of vy0 using step_off in x direction from Li

    Note: Typical value of pert_x in Earth-Moon CR3BP is 0.01 [nd]

    Dhruv Jain, Feb 25 2022

    Parameters
    ----------
    mu :  float, M2/(M1+M2)
        M1 and M2 are mass of Primary Bodies and M2<M1
    pert_x : float
        Step off from Li in x direction, CR3BP Barycentered frame
    lib_num : int, optional
        1: L1, 2: L2, 3: L3. The default is 1.

    Returns
    -------
    initial guess: numpy ndarray (6x1)
        Initial guess of a lyapunov orbit about a collinear libration point calculated
        using CR3BP in-plane linear variational equation
    """

    # Check if lib_num is valid
    if lib_num not in [1, 2, 3]:
        print("Function only works for collinear libration points: 1, 2 or 3")
        return 0

    # Calculate the 5 libration points
    lib_loc = lib_pt_loc(mu)
    li = lib_loc[lib_num - 1, :]  # 0 for L1, 1 for L2, 2 for L3
    li_state = np.array([li[0], 0, 0, 0, 0, 0])

    d, r = rel_dist_cr3bp(mu, li_state)
    Uxx, Uyy, Uzz, Uxy, Uxz, Uyz = uii_partials_cr3bp(
        mu, li_state, d, r
    )  # Hetian of  Pesudo-potenital of CR3BP

    # Compute roots of characteristic polynomial of In-plane Linear Variational Equation of CR3BP
    # Characteristic polynomial: Lam^4 + (4 -Uxx - Uyy)*Lam^2 + Uxx*Uyy = 0
    # Characteristic polynomial: lam^2 + 2*b1*lam - b2^2 = 0
    b1 = 2 - (Uxx + Uyy) / 2
    b2_2 = -Uxx * Uyy  # b2^2

    s = (b1 + (b1**2 + b2_2) ** 0.5) ** 0.5
    b3 = (s**2 + Uxx) / (2 * s)

    # Step_off from Li to compute linear approximation of vy0
    enz0i = np.array([pert_x, 0, 0])
    enzdot_0i = np.array([enz0i[1] * s / b3, -b3 * enz0i[0] * s, 0])
    delta_x0 = np.concatenate((enz0i, enzdot_0i))

    initial_guess = li_state + delta_x0

    return initial_guess
