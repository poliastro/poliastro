"""
@author: Dhruv Jain, Multi-Body Dynaminitial_guesss Research Group, MSAAE Purdue University

Objectve: Calculates the position [nd] of 5 libration points of a CR3BP system
"""
import numpy as np
from astropy import units as u


def lib_pt_loc(SysChars, conv_tol=1e-12):
    """Computes libration points position [nd] for a CR3BP system

    Parameters
    ----------
    SysChars: object
        Object of Class SystemChars
    conv_tol: float
        convergence tolerance for Newton-Raphson Method

    Returns
    -------
    lib_loc: numpy ndarray (5x3)
        5 Libration Points, [nd]
    """

    mu = SysChars.mu

    lib_loc = np.zeros((5, 3))

    # 5th degree polynomial of L1, L2 and L3
    f_lib123_coeffs = np.array(
        [
            [1, mu - 3, 3 - 2 * mu, -mu, 2 * mu, -mu],
            [1, 3 - mu, 3 - 2 * mu, -mu, -2 * mu, -mu],
            [1, 2 + mu, 1 + 2 * mu, mu - 1, 2 * mu - 2, -1 + mu],
        ]
    )

    # First-order derivative of the polyomial defined in f_lib123 with respect to 'Li,x'; e.g. L1,x ; L2,x ; L3,x
    df_lib123_coeffs = np.array(
        [
            [0, 5, 4 * (mu - 3), 3 * (3 - 2 * mu), 2 * -mu, 2 * mu],
            [0, 5, 4 * (3 - mu), 3 * (3 - 2 * mu), 2 * -mu, -2 * mu],
            [0, 5, 4 * (2 + mu), 3 * (1 + 2 * mu), 2 * (mu - 1), 2 * mu - 2],
        ]
    )

    initial_guess123 = np.array([0.9, 1.1, -1])  # Initial guess for L1, L2, L3

    # computes L1
    lib_loc[0, 0] = (
        1
        - mu
        - newton_raphson_lib_calc(
            initial_guess123[0],
            f_lib123_coeffs[0, :],
            df_lib123_coeffs[0, :],
            conv_tol,
        )
    )

    # computes L2
    lib_loc[1, 0] = (
        1
        - mu
        + newton_raphson_lib_calc(
            initial_guess123[1],
            f_lib123_coeffs[1, :],
            df_lib123_coeffs[1, :],
            conv_tol,
        )
    )

    # computes L3
    lib_loc[2, 0] = -mu - newton_raphson_lib_calc(
        initial_guess123[2],
        f_lib123_coeffs[2, :],
        df_lib123_coeffs[2, :],
        conv_tol,
    )

    # L4, analytical solution
    lib_loc[3, :] = [
        0.5 - mu,
        3**0.5 / 2,
        0,
    ]

    # L5, analytical solution
    lib_loc[4, :] = [
        0.5 - mu,
        -(3**0.5) / 2,
        0,
    ]

    # # Define custome units, [nd]: non-dimensional unit for distance
    L_ND = u.def_unit("dist_nd", SysChars.lstar)
    # # u.add_enabled_units(L_ND)

    return lib_loc * L_ND


def newton_raphson_lib_calc(
    initial_guess, func_coeffs, dfunc_coeffs, conv_tol
):
    """Uses Newton-Raphson Method to compute a zero of a function

    Parameters
    ----------
    initial_guess: float
        Approximate value of a zero of the function 'func'
    func_coeffs: numpy ndarray (6x1)
        Coeffecients of a polynomial function whose zero is to be computed
        [C_n, C_n-1, C_n-2, ...., C_0]
        Setup: C_n * x^n + C_n-1 * x^(n-1) ..... C0 * x^0
    dfunc_coeffs: numpy ndarray (6x1)
        Coeffecients of the derivative of function 'func' w.r.t 'x' whose zero is to be computed
        [C_n, C_n-1, C_n-2, ...., C_0]
        Setup: C_n * x^(n-1) + C_n-1 * x^(n-2) ..... C0 * x^0
    conv_tol: float
        convergence tolerance for Newton-Raphson Method

    Returns
    -------
    initial_guess: float
        zero of function 'func_coeffs' within conv_tol
    """

    initial_guess = np.array(
        [initial_guess]
    )  # To increase dimension to 1, which is the minimum required dimension for np.vander input

    val = np.vander(initial_guess, N=6)
    h = np.dot(val, func_coeffs) / np.dot(val, dfunc_coeffs)

    # Using Newton-Raphson Method
    while abs(h) >= conv_tol:

        val = np.vander(initial_guess, N=6)
        h = np.dot(val, func_coeffs) / np.dot(val, dfunc_coeffs)
        initial_guess = initial_guess - h

    return initial_guess
