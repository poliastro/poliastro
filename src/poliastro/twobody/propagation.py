"""Propagation algorithms

"""
import functools

import numpy as np
import scipy as sp
from astropy import time, units as u
from astropy.coordinates import CartesianDifferential, CartesianRepresentation
from scipy.integrate import solve_ivp

from poliastro.core.propagation import (
    func_twobody,
    kepler as kepler_fast,
    mean_motion as mean_motion_fast,
)
from poliastro.core.util import norm
from poliastro.integrators import DOP835


def cowell(k, r, v, tofs, rtol=1e-11, *, ad=None, **ad_kwargs):
    """Propagates orbit using Cowell's formulation.

    Parameters
    ----------
    k : ~astropy.units.Quantity
        Standard gravitational parameter of the attractor.
    r : ~astropy.units.Quantity
        Position vector.
    v : ~astropy.units.Quantity
        Velocity vector.
    tofs : ~astropy.units.Quantity
        Array of times to propagate.
    rtol : float, optional
        Maximum relative error permitted, default to 1e-10.
    ad : function(t0, u, k), optional
         Non Keplerian acceleration (km/s2), default to None.

    Returns
    -------
    rr : ~astropy.units.Quantity
        Propagated position vectors.
    vv : ~astropy.units.Quantity
        Propagated velocity vectors.

    Raises
    ------
    RuntimeError
        If the algorithm didn't converge.

    Note
    -----
    This method uses a Dormand & Prince method of order 8(5,3) available
    in the :py:class:`poliastro.integrators` module. If multiple tofs
    are provided, the method propagates to the maximum value and
    calculates the other values via dense output

    """
    k = k.to(u.km ** 3 / u.s ** 2).value
    x, y, z = r.to(u.km).value
    vx, vy, vz = v.to(u.km / u.s).value
    tofs = tofs.to(u.s).value

    u0 = np.array([x, y, z, vx, vy, vz])

    # Set the non Keplerian acceleration
    if ad is None:

        def ad(t0, u_, k_):
            return 0, 0, 0

    f_with_ad = functools.partial(func_twobody, k=k, ad=ad, ad_kwargs=ad_kwargs)

    result = solve_ivp(
        f_with_ad,
        (0, max(tofs)),
        u0,
        rtol=rtol,
        atol=1e-12,
        method=DOP835,
        dense_output=True,
    )
    if not result.success:
        raise RuntimeError("Integration failed")

    rrs = []
    vvs = []
    for i in range(len(tofs)):
        y = result.sol(tofs[i])
        rrs.append(y[:3])
        vvs.append(y[3:])

    return rrs * u.km, vvs * u.km / u.s


def mean_motion(k, r, v, tofs, **kwargs):
    """Propagates orbit using Cowell's formulation.

    Parameters
    ----------
    k : ~astropy.units.Quantity
        Standard gravitational parameter of the attractor.
    r : ~astropy.units.Quantity
        Position vector.
    v : ~astropy.units.Quantity
        Velocity vector.
    tofs : ~astropy.units.Quantity
        Array of times to propagate.

    Returns
    -------
    rr : ~astropy.units.Quantity
        Propagated position vectors.
    vv : ~astropy.units.Quantity
        Propagated velocity vectors.

    """
    k = k.to(u.km ** 3 / u.s ** 2).value
    r0 = r.to(u.km).value
    v0 = v.to(u.km / u.s).value
    tofs = tofs.to(u.s).value

    results = [mean_motion_fast(k, r0, v0, tof) for tof in tofs]
    # TODO: Rewrite to avoid iterating twice
    return (
        [result[0] for result in results] * u.km,
        [result[1] for result in results] * u.km / u.s,
    )


def kepler(k, r, v, tofs, numiter=350, **kwargs):
    """Propagates Keplerian orbit.

    Parameters
    ----------
    k : ~astropy.units.Quantity
        Standard gravitational parameter of the attractor.
    r : ~astropy.units.Quantity
        Position vector.
    v : ~astropy.units.Quantity
        Velocity vector.
    tofs : ~astropy.units.Quantity
        Array of times to propagate.
    numiter : int, optional
        Maximum number of iterations, default to 35.

    Returns
    -------
    rr : ~astropy.units.Quantity
        Propagated position vectors.
    vv : ~astropy.units.Quantity
        Propagated velocity vectors.

    Raises
    ------
    RuntimeError
        If the algorithm didn't converge.

    Note
    -----
    This algorithm is based on Vallado implementation, and does basic Newton
    iteration on the Kepler equation written using universal variables. Battin
    claims his algorithm uses the same amount of memory but is between 40 %
    and 85 % faster.

    """
    k = k.to(u.km ** 3 / u.s ** 2).value
    r0 = r.to(u.km).value
    v0 = v.to(u.km / u.s).value
    tofs = tofs.to(u.s).value

    results = [_kepler(k, r0, v0, tof, numiter=numiter) for tof in tofs]
    # TODO: Rewrite to avoid iterating twice
    return (
        [result[0] for result in results] * u.km,
        [result[1] for result in results] * u.km / u.s,
    )


def _kepler(k, r0, v0, tof, *, numiter):
    # Compute Lagrange coefficients
    f, g, fdot, gdot = kepler_fast(k, r0, v0, tof, numiter)

    assert np.abs(f * gdot - fdot * g - 1) < 1e-5  # Fixed tolerance

    # Return position and velocity vectors
    r = f * r0 + g * v0
    v = fdot * r0 + gdot * v0

    return r, v


def _compute_C(A):
    """ This function returns C given A. """

    # Taken directly from Battin 1999
    c_coeff_frac = [
        8 / 175,
        8 / 525,
        1896 / 336875,
        28744 / 13138125,
        26248 / 29859375,
        134584 / 372246875,
        129802986344 / 857740555546875,
        55082676858 / 857740555546875,
        687061097149992 / 24934041427216796875,
        8033038585237352 / 673219118534853515625,
        2892031498456202296 / 555405772791254150390625,
        4093458329892811912 / 1789640823438485595703125,
        10507274811793509548037032 / 10398126371321703046014404296875,
        25259289756593760965336 / 56299726292037542523193359375,
        329740626632896883140474984552 / 1648016378801395585267899627685546875,
        3277387142806000848933442276712 / 36585963609390981992947371734619140625,
        70279162673586836547221944640968 / 1746148263175478686027033650970458984375,
        2556121670455925738199518346408 / 140904073482222380998066233062744140625,
    ]

    # Taken directly from: UNIFIED SYMBOLIC ALGORITHM OF GAUSS METHOD FOR NEAR-PARABOLIC ORBITS
    c_coeff_float = [
        4.571428571428571e-2,
        1.523809523809524e-2,
        5.628200371057514e-3,
        2.187831216402645e-3,
        8.790538984824699e-4,
        3.615450096122365e-4,
        1.513312918511131e-4,
        6.421834259764084e-5,
        2.755514380432645e-5,
        1.193227934869094e-5,
        5.207060567487393e-6,
        2.287307193869180e-6,
        1.010496933444937e-6,
        4.486574166554373e-7,
        2.000833431477894e-7,
        8.958045161245261e-8,
        4.024810730892910e-8,
        1.814086425811123e-8,
    ]

    C = 0
    for i, c_i in enumerate(c_coeff_float, 2):
        C += c_i * A ** i

    return C


def _gauss(k, r0, v0, tof, numiter=350, rtol=1e-12, **kwargs):
    """ This is improved Gauss algorithm for solving near parabolic orbits.
    It was developed by Battin and Fill after publishing their paper 'Extension
    of Gauss method for the solution of Keplers Equation'.

    Parameters
    ----------
    k: astropy.quantities.Quantity
        Gravitational parameter
    r0: np.array
        Initial position vector
    v0: np.array
        Initial velocity vector
    tof: astropy.quantity.Quantity
        Time of flight to propagate
    numiter: int
        Number of iterations for iterative method
    rtol: float
        Tolerance for the method

    Returns
    -------
    r: np.array
        Final position vector
    v: np.aray
        Final velocity vector

    Note
    ----
    This algorithm was directly taken from Battin 1999, "An Introduction to
    the Mathematics and Methods of Astrodynamics"

    """

    # Just storing this values for further computation
    state_vector = np.array([[r0[0], r0[1], r0[2]], [v0[0], v0[1], v0[2]]])

    # Compute dot product r · v
    dot_rv = np.dot(r0, v0)

    # Working with modulus of r0 and v0
    r0 = norm(r0)
    v0 = norm(v0)

    # Initial parameters for the method
    delta0 = r0 * v0 ** 2 / k
    beta0 = np.sqrt(k / r0 ** 3)
    phi0 = dot_rv / np.sqrt(k * r0)
    T_prime = 3 * beta0 * tof

    # Main variables for the numerical method
    A = 0
    A_m = 0.3
    n = 0

    # Iteration begins
    while n < numiter:

        # Computing main variables of algorithm
        # These ones need to be inside of the loop, since its
        # value is modified if A >= A_m
        phi0_prime = 0.5 * phi0
        phi0_2prime = 0.5 * phi0_prime
        gamma0 = 4 * (1 / 2 - 1 / 4 * delta0)
        gamma0_prime = 1 / 4 * gamma0
        gamma0_2prime = 1 - gamma0
        ji0 = 2 * (9 / 40 * delta0 - 1 / 5)
        ji0_prime = 1 / 2 * ji0

        # Algorithm begins iterations
        C = _compute_C(A)

        # Equation (5.93)
        alpha = 1 + 1 / 5 * A + C
        alpha_prime = 1 - 3 / 5 * A

        # Equation (5.95)
        psi = alpha * T_prime
        eta = alpha_prime * alpha_prime * psi / (alpha - A)

        # Equation (5.96)
        epsilon = 1 + eta * phi0_prime
        b = np.abs(epsilon + eta * (phi0_2prime + ji0_prime * psi))

        # TODO: Use custom Newton-Raphson method
        def f(x, epsilon, b):
            # Equation (5.83)
            return x ** 3 - 3 * epsilon * x - 2 * b

        x = sp.optimize.newton(f, b, args=(epsilon, b), tol=1e-12)

        # Equation (5.97)
        if (epsilon + eta * (phi0_2prime + ji0_prime * psi)) >= 0:
            xi = 1 + x
        else:
            xi = 1 - x

        theta = psi * eta / (xi ** 2)

        # Equation (5.97)
        A_new = gamma0_prime * theta

        # Each time A_new is computed a test to check if is lower than A_m
        # If it is higher, a recomputation of the initial parameters for the
        # algorithm is done.
        if np.abs(A_new) >= A_m:

            A_old = A_new

            # We set the magnitude of A ot the value of A_m
            A_new = np.sign(A_new) * A_m
            C = _compute_C(A_new)
            alpha_m = 1 + 1 / 5 * A_new + C
            alpha_prime_m = 1 - 3 / 5 * A_new

            # Equation (5.101)
            theta_m = A_new / gamma0_prime
            lamda_m = theta_m / 2 / alpha_m
            B = alpha_m * alpha_prime_m / np.sqrt(alpha_m - A_new)
            Tm_prime = B * np.sqrt(theta_m) * (3 + ji0 * theta_m) + 3 * phi0 * lamda_m

            # Equation (5.98) for maximum case
            lamda_m = theta_m / 2 / alpha_m
            k_m = alpha_prime_m * Tm_prime / xi

            # Equation (5.99)
            rho_m = 1 + gamma0_2prime * lamda_m + phi0 * k_m

            # Equation (5.100)
            matrix_phi = np.array(
                [
                    [1 - lamda_m, (k_m + phi0 * lamda_m) / beta0],
                    [-beta0 * k_m / rho_m, 1 - lamda_m / rho_m],
                ]
            )

            state_vector = np.matmul(matrix_phi, state_vector)

            # Equation (5.102)
            r_m = rho_m * r0
            phi_m = (phi0 * (1 - gamma0 * lamda_m) + gamma0_2prime * k_m) / np.sqrt(
                rho_m
            )
            beta_m = beta0 / rho_m / np.sqrt(rho_m)
            delta_m = 2 - rho_m * gamma0

            # We replace now the parameters of the algorithm
            T_prime = T_prime - Tm_prime
            r0 = r_m
            beta0 = beta_m
            phi0 = phi_m
            delta0 = delta_m

            # Increase one iteration
            A = A_old
            n += 1

        elif np.abs(A - A_new) <= rtol:
            # Solution has converged and we can return the result

            # Equation (5.98)
            lamda_ = theta / 2 / alpha
            k = alpha_prime * T_prime / xi

            # Equation (5.99)
            rho = 1 + gamma0_2prime * lamda_ + phi0 * k

            # Equation (5.100)
            matrix_phi = np.array(
                [
                    [1 - lamda_, (k + phi0 * lamda_) / beta0],
                    [-beta0 * k / rho, 1 - lamda_ / rho],
                ]
            )

            state_vector = np.matmul(matrix_phi, state_vector)

            return state_vector
        else:
            # Increase one iteration
            A = A_new
            n += 1


def gauss(k, r0, v0, time_of_flight, **kwargs):

    k = k.to(u.m ** 3 / u.s ** 2).value
    r0 = r0.to(u.m).value
    v0 = v0.to(u.m / u.s).value
    time_of_flight = time_of_flight.to(u.s).value

    results = [_gauss(k, r0, v0, tof, **kwargs) for tof in time_of_flight]

    return (
        [result[0] for result in results] * u.m,
        [result[1] for result in results] * u.m / u.s,
    )


def propagate(orbit, time_of_flight, *, method=mean_motion, rtol=1e-10, **kwargs):
    """Propagate an orbit some time and return the result.

    Parameters
    ----------
    orbit : ~poliastro.twobody.Orbit
        Orbit object to propagate.
    time_of_flight : ~astropy.time.TimeDelta
        Time of propagation.
    method : callable, optional
        Propagation method, default to mean_motion.
    rtol : float, optional
        Relative tolerance, default to 1e-10.

    """
    # Use the highest precision we can afford
    # np.atleast_1d does not work directly on TimeDelta objects
    jd1 = np.atleast_1d(time_of_flight.jd1)
    jd2 = np.atleast_1d(time_of_flight.jd2)
    time_of_flight = time.TimeDelta(jd1, jd2, format="jd", scale=time_of_flight.scale)

    rr, vv = method(
        orbit.attractor.k, orbit.r, orbit.v, time_of_flight.to(u.s), rtol=rtol, **kwargs
    )

    # TODO: Turn these into unit tests
    assert rr.ndim == 2
    assert vv.ndim == 2

    cartesian = CartesianRepresentation(
        rr, differentials=CartesianDifferential(vv, xyz_axis=1), xyz_axis=1
    )

    return cartesian
