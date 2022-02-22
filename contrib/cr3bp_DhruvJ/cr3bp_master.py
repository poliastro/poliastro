"""
Created on Jan 8 21:42:49 2022

@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University
        dhruvj9922@gmail.com
        
Objective: This file contains key functions required to numerically integrate a 
    state defined about the barycenter of two primaries in the Circular Restricted
    Three Body Problem (CR3BP) model
    
    Features: 
        1. Integrate CR3BP EOMs
        2. Integrate CR3BP EOMs + Compute the State Transition Matrix for each state
        3. Optional events function is added to numerical integrator to track when states move from -y to +y region and/or visa versa
        4. Computes accelration, first-derivate of pseudo-potential, and second-derivative of pseudo-potenital terms terms

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
from scipy.integrate import solve_ivp


def prop_cr3bp(mu, ic, tf, tol=1e-12, stm_bool=0, xcross_cond=0, int_method="DOP853"):
    """Numerically Integrate Circular Restricted Three-Body Problem EOMs
    Dhruv Jain, Jan 8 2022

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
    tol : float, optional
        Absolute = Relative Integration Tolerance
        The default is 1e-12.
    stm_bool : boolean, optional
        0: CR3BP EOM Integration
        1: CR3BP + STM EOM Integration
        The default is 0.
    xcross_cond : int, optional
        0 => No, else x crossing condition

        DESCRIPTION. The default is 0.
    int_method : string, optional
        Specify integration scheme: 'DOP853' or 'LSODA'
        The default is 'DOP853'.

    Returns
    -------
     results : Dictionary
         't': time history, [nd]
         'states': state history, [:,i] => 6 states @ t = t_i, [nd]
         'yevents': states at prescribed event, [nd]
         'tevents': time stamp of yevents, [nd]
         'stm': STM element history, [:,:,i] => STM @ t = t_i
    """

    # Check to see if any state of I.C. is complex and then accordingly set the datatype
    # This is done to ease the implementation of Complex step derivative formulation to compute numerical partials in the future
    if all(np.imag(imag_check) == 0 for imag_check in ic):
        datatype = np.float64
        ic = np.real(ic)
    else:
        datatype = np.complex
        if int_method == "LSODA":
            print(int_method + " cannot integrate in complex domain")
            return 0

    # Check integration scheme
    if int_method != "DOP853" and int_method != "LSODA":
        print(
            "Please DROP853 or LSODA as the integration schemes, other schemes may not be compatible with the setup"
        )
        return 0

    # Accept 6 states and append 36 inital states if stm_bool = 1
    # Runs even if all 42 sates are given
    if len(ic) == 6 and stm_bool == 1:
        ic = np.concatenate(
            (ic, np.identity(6).flatten())
        )  # Appends the IC for STM: I_6x6
        print(ic)
    elif len(ic) == 42:
        stm_bool = 1  # To make sure stm_boolis set to 1 if all 42 states are passed
    if len(ic) == 6:
        stm_bool == 0
    elif len(ic) != 6 or len(ic) != 42:
        print(
            "Initial conditions are niether of length 6 nor 42, recheck the input, len is"
            + str(len(ic))
        )
        return 0

    # Skeleton function for CR3BP + STM Numerical Integration
    def NDDE_STM(t, x):

        n = 1
        d, r = rel_dist_cr3bp(mu, x)  # Relative position vectors

        dx = np.empty((len(x),), dtype=datatype)

        # CR3BP: 6 states

        # CR3BP: 6 states
        dx[0] = x[3]
        dx[1] = x[4]
        dx[2] = x[5]
        dx[3] = (
            2 * n * x[4]
            + n**2 * x[0]
            - (1 - mu) * (x[0] + mu) / d**3
            - mu * (x[0] - 1 + mu) / r**3
        )
        dx[4] = (
            -2 * n * x[3]
            + n**2 * x[1]
            - (1 - mu) * x[1] / d**3
            - mu * x[1] / r**3
        )
        dx[5] = -(1 - mu) * x[2] / d**3 - mu * x[2] / r**3

        # STM: 36 States
        if len(x) == 42:
            Uxx, Uyy, Uzz, Uxy, Uxz, Uyz = uii_partials_cr3bp(
                mu, x, d, r
            )  # Hetian of  Pesudo-potenital of CR3BP
            At = np.array(
                [
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1],
                    [Uxx, Uxy, Uxz, 0, 2, 0],
                    [Uxy, Uyy, Uyz, -2, 0, 0],
                    [Uxz, Uyz, Uzz, 0, 0, 0],
                ]
            )

            phi = np.zeros((6, 6), dtype=datatype)

            # Assingment of x[] for phi, STM
            count = 6
            for i in range(6):
                for j in range(6):
                    phi[i][j] = x[count]
                    count = count + 1

            phi_d = np.matmul(At, phi)  # STM_dot = A(t)*STM

            dx[6:42] = phi_d.flatten()

        return dx

    t0 = 0  # By default set initial integration time

    def xcross(t, y):
        """Track y position state for events function during integration"""
        return y[1]

    if (
        xcross_cond == 1
    ):  # Track events when crossising -y to +y region or +y to -y region, that is crossing XZ or XY plane
        xcross.terminal = False
        xcross.direction = 0
        fun = solve_ivp(
            NDDE_STM, [t0, tf], ic, method=int_method, events=xcross, rtol=tol, atol=tol
        )
    elif xcross_cond == 2:  # Track events when crossing -y to +y region
        xcross.terminal = True
        xcross.direction = 1
        fun = solve_ivp(
            NDDE_STM, [t0, tf], ic, method=int_method, events=xcross, rtol=tol, atol=tol
        )
    else:  # Does not track any events
        fun = solve_ivp(NDDE_STM, [t0, tf], ic, method=int_method, rtol=tol, atol=tol)

    # Save data to dictionary
    results = save_prop_data_cr3p(fun, stm_bool, datatype)

    return results


def save_prop_data_cr3p(fun, stm_bool, datatype):
    """Saves Numerical Integration results to a dictionary
    Dhruv Jain, Jan 9 2022

    Parameters
    ----------
    fun : Object
        solve_ivp object
    stm_bool : boolean, optional
        0: CR3BP EOM Integration
        1: CR3BP + STM EOM Integration
    datatype : datatype
        Can be np.float64 or np.complex

    Returns
    -------
    results : Dictionary
        't': time history, [nd]
        'states': state history, [:,i] => 6 states @ t = t_i, [nd]
        'yevents': states at prescribed event, [nd]
        'tevents': time stamp of yevents, [nd]
        'stm': STM element history, [:,:,i] => STM @ t = t_i

    """
    # Save the 6 states
    t = fun.t
    x = fun.y[0]
    y = fun.y[1]
    z = fun.y[2]
    vx = fun.y[3]
    vy = fun.y[4]
    vz = fun.y[5]

    # Save events function results
    tevents = fun.t_events
    yevents = fun.y_events

    # Save STM elements
    if stm_bool == 1:
        count = 0
        stm_vals = np.zeros((6, 6, len(t)), dtype=datatype)
        for t1 in range(6):
            for t2 in range(6):
                for t3 in range(len(t)):
                    stm_vals[t1, t2, t3] = fun.y[count + 6][t3]
                count = count + 1

    results = {}
    results["t"] = np.asarray(t)
    results["states"] = np.concatenate(
        [
            [np.asarray(x)],
            [np.asarray(y)],
            [np.asarray(z)],
            [np.asarray(vx)],
            [np.asarray(vy)],
            [np.asarray(vz)],
        ]
    ).T
    if stm_bool == 1:
        results["stm"] = stm_vals
    results["tevents"] = tevents
    results["yevents"] = yevents

    return results


def rel_dist_cr3bp(mu, states):
    """Compute distance between a satellite(P3) defined in P1-P2 barycenter
        to P1 and P2 in CR3BP
    Dhruv Jain, Jan 9 2022

    Parameters
    ----------
    mu :  float, M2/(M1+M2)
        M1 and M2 are mass of Primary Bodies and M2<M1
    states : numpy ndarray (6x1)
        States are defined about the barycenter of the two primaries, P1 and P2

    Returns
    -------
    d : float64/complex128
        distance between P3 and P1
    r : float64/complex128
        distance between P3 and P2
    """
    d = ((states[0] + mu) ** 2 + states[1] ** 2 + states[2] ** 2) ** 0.5
    r = ((states[0] - 1 + mu) ** 2 + states[1] ** 2 + states[2] ** 2) ** 0.5

    return d, r


def uii_partials_cr3bp(mu, states, d, r):
    """Compute second-derivate of the pseudo-potenital of the CR3BP EOMs
    Dhruv Jain, Jan 10 2022

    Parameters
    ----------
    mu :  float, M2/(M1+M2)
        M1 and M2 are mass of Primary Bodies and M2<M1
    states : numpy ndarray (6x1)
        States are defined about the barycenter of the two primaries, P1 and P2
    d : float64/complex128
        distance between P3 and P1
    r : float64/complex128
        distance between P3 and P2

    Returns
    -------
    Uxx, Uyy, Uzz, Uxy, Uxz, Uyz: float64/complex128
        Second-derivaitves of the pseudo-potenial of CR3BP

    """
    # Second order partials of U for A(t) matrix (6x6)
    Uxx = (
        1
        - (1 - mu) / d**3
        - mu / r**3
        + 3 * (1 - mu) * (states[0] + mu) ** 2 / d**5
        + 3 * mu * (states[0] - 1 + mu) ** 2 / r**5
    )
    Uyy = (
        1
        - (1 - mu) / d**3
        - mu / r**3
        + 3 * (1 - mu) * states[1] ** 2 / d**5
        + 3 * mu * states[1] ** 2 / r**5
    )
    Uzz = (
        -(1 - mu) / d**3
        - mu / r**3
        + 3 * (1 - mu) * states[2] ** 2 / d**5
        + 3 * mu * states[2] ** 2 / r**5
    )
    Uxy = (
        3 * (1 - mu) * (states[0] + mu) * states[1] / d**5
        + 3 * mu * (states[0] - 1 + mu) * states[1] / r**5
    )
    Uxz = (
        3 * (1 - mu) * (states[0] + mu) * states[2] / d**5
        + 3 * mu * (states[0] - 1 + mu) * states[2] / r**5
    )
    Uyz = (
        3 * (1 - mu) * states[1] * states[2] / d**5
        + 3 * mu * states[1] * states[2] / r**5
    )

    return Uxx, Uyy, Uzz, Uxy, Uxz, Uyz


def ui_partials_acc_cr3bp(mu, states):
    """Compute first-derivateive of pseudo-potenital terms of the CR3BP EOMs and acceleration terms
    Dhruv Jain, Jan 10 2022

    Parameters
    ----------
    mu :  float, M2/(M1+M2)
        M1 and M2 are mass of Primary Bodies and M2<M1
    states : numpy ndarray (6x1)
        States are defined about the barycenter of the two primaries, P1 and P2
    """
    d, r = rel_dist_cr3bp(mu, states)
    Ux = (
        states[0]
        - (1 - mu) / d**3 * (states[0] + mu)
        - mu / r**3 * (states[0] - 1 + mu)
    )
    Uy = states[1] - (1 - mu) / d**3 * states[1] - mu / r**3 * states[1]
    Uz = -(1 - mu) / d**3 * states[2] - mu / r**3 * states[2]
    ax = 2 * states[4] + Ux
    ay = -2 * states[3] + Uy
    az = Uz

    return Ux, Uy, Uz, ax, ay, az
