# -*- coding: utf-8 -*-
"""
Created on 7 Jan 22 12:24:03 
Updated 20 Mar 2022

@author: Dhruv Jain, Multi-Body Dynaminitial_guesss Research Group, MSAAE Purdue University
        dhruvj9922@gmail.com

Objectve: Calculate the position [nd] of 5 libration points for a system in CR3BP
"""
import numpy as np


def lib_pt_loc(sys_chars_vals, tolerance=1e-12):
    """Computes Non-Dimensionalized Libration Points Location for P1-P2 system
    Parameters
    ----------
    sys_chars_vals: object
        Object of Class sys_char

    tolerance: float
        convergence tolerance for Newton-Raphson Method

    Returns
    -------
    lib_loc: numpy ndarray (5x3)
        5 Libration Points, [nd]
    """
    mu = sys_chars_vals.mu

    lib_loc = np.zeros((5, 3))
    lib_loc[3, :] = [
        0.5 - mu,
        3**0.5 / 2,
        0,
    ]  # L4, analytical_guessal solution known
    lib_loc[4, :] = [
        0.5 - mu,
        -(3**0.5) / 2,
        0,
    ]  # L5, analytical solution known

    # 5th degree polynomial of L1, L2 and L3
    f_lib = np.array(
        [
            [1, mu - 3, 3 - 2 * mu, -mu, 2 * mu, -mu],
            [1, 3 - mu, 3 - 2 * mu, -mu, -2 * mu, -mu],
            [1, 2 + mu, 1 + 2 * mu, mu - 1, 2 * mu - 2, -1 + mu],
        ]
    )

    # First-order derivative of the polyomial defined in f_lib
    fd_lib = np.array(
        [
            [0, 5, 4 * (mu - 3), 3 * (3 - 2 * mu), 2 * -mu, 2 * mu],
            [0, 5, 4 * (3 - mu), 3 * (3 - 2 * mu), 2 * -mu, -2 * mu],
            [0, 5, 4 * (2 + mu), 3 * (1 + 2 * mu), 2 * (mu - 1), 2 * mu - 2],
        ]
    )

    initial_guess = np.array([0.9, 1.1, -1])

    for i in range(3):
        val = np.vander([initial_guess[i]], 6)
        h = np.dot(val, f_lib[i, :]) / np.dot(val, fd_lib[i, :])
        while abs(h) >= tolerance:
            val = np.vander([initial_guess[i]], 6)
            h = np.dot(val, f_lib[i, :]) / np.dot(val, fd_lib[i, :])
            lib_loc[i, 0] = initial_guess[i] - h

            initial_guess[i] = lib_loc[i, 0]

        if i == 0:
            lib_loc[i, 0] = 1 - mu - lib_loc[i, 0]
        elif i == 1:
            lib_loc[i, 0] = 1 - mu + lib_loc[i, 0]
        elif i == 2:
            lib_loc[i, 0] = -mu - lib_loc[i, 0]

    return lib_loc
