# coding: utf-8
"""Propagation algorithms.

"""
import numpy as np

from scipy.integrate import ode

from astropy import units as u

from poliastro.jit import jit
from poliastro.stumpff import c2, c3


def func_twobody(t0, u_, k, ad):
    """Differential equation for the initial value two body problem.

    This function follows Cowell's formulation.

    Parameters
    ----------
    t0 : float
        Time.
    u_ : ndarray
        Six component state vector [x, y, z, vx, vy, vz] (km, km/s).
    k : float
        Standard gravitational parameter.
    ad : function(t0, u, k)
         Non Keplerian acceleration (km/s2).

    """
    ax, ay, az = ad(t0, u_, k)

    x, y, z, vx, vy, vz = u_
    r3 = (x**2 + y**2 + z**2)**1.5

    du = np.array([
        vx,
        vy,
        vz,
        -k * x / r3 + ax,
        -k * y / r3 + ay,
        -k * z / r3 + az
    ])
    return du


def cowell(k, r0, v0, tof, rtol=1e-10, *, ad=None, callback=None, nsteps=1000):
    """Propagates orbit using Cowell's formulation.

    Parameters
    ----------
    k : float
        Gravitational constant of main attractor (km^3 / s^2).
    r0 : array
        Initial position (km).
    v0 : array
        Initial velocity (km).
    ad : function(t0, u, k), optional
         Non Keplerian acceleration (km/s2), default to None.
    tof : float
        Time of flight (s).
    rtol : float, optional
        Maximum relative error permitted, default to 1e-10.
    nsteps : int, optional
        Maximum number of internal steps, default to 1000.
    callback : callable, optional
        Function called at each internal integrator step.

    Raises
    ------
    RuntimeError
        If the algorithm didn't converge.

    Notes
    -----
    This method uses a Dormand & Prince method of order 8(5,3) available
    in the ``scipy.integrate.ode`` module.

    """
    x, y, z = r0
    vx, vy, vz = v0
    u0 = np.array([x, y, z, vx, vy, vz])

    # Set the non Keplerian acceleration
    if ad is None:
        ad = lambda t0, u_, k_: (0, 0, 0)

    # Set the integrator
    rr = ode(func_twobody).set_integrator('dop853', rtol=rtol, nsteps=nsteps)
    rr.set_initial_value(u0)  # Initial time equal to 0.0
    rr.set_f_params(k, ad)  # Parameters of the integration
    if callback:
        rr.set_solout(callback)

    # Make integration step
    rr.integrate(tof)

    if rr.successful():
        r, v = rr.y[:3], rr.y[3:]
    else:
        raise RuntimeError("Integration failed")

    return r, v


def kepler(k, r0, v0, tof, rtol=1e-10, *, numiter=35):
    """Propagates Keplerian orbit.

    Parameters
    ----------
    k : float
        Gravitational constant of main attractor (km^3 / s^2).
    r0 : array
        Initial position (km).
    v0 : array
        Initial velocity (km).
    tof : float
        Time of flight (s).
    rtol : float, optional
        Maximum relative error permitted, default to 1e-10.
    numiter : int, optional
        Maximum number of iterations, default to 35.

    Raises
    ------
    RuntimeError
        If the algorithm didn't converge.

    Notes
    -----
    This algorithm is based on Vallado implementation, and does basic Newton
    iteration on the Kepler equation written using universal variables. Battin
    claims his algorithm uses the same amount of memory but is between 40 %
    and 85 % faster.

    """
    # Compute Lagrange coefficients
    f, g, fdot, gdot = _kepler(k, r0, v0, tof, numiter, rtol)

    assert np.abs(f * gdot - fdot * g - 1) < 1e-5  # Fixed tolerance

    # Return position and velocity vectors
    r = f * r0 + g * v0
    v = fdot * r0 + gdot * v0

    return r, v


def propagate(orbit, time_of_flight, *, method=kepler, rtol=1e-10, **kwargs):
    """Propagate an orbit some time and return the result.

    """
    r, v = method(orbit.attractor.k.to(u.km ** 3 / u.s ** 2).value,
                  orbit.r.to(u.km).value, orbit.v.to(u.km / u.s).value,
                  time_of_flight.to(u.s).value,
                  rtol=rtol,
                  **kwargs)
    return orbit.from_vectors(orbit.attractor, r * u.km, v * u.km / u.s, orbit.epoch + time_of_flight)


@jit
def _kepler(k, r0, v0, tof, numiter, rtol):
    # Cache some results
    dot_r0v0 = np.dot(r0, v0)
    norm_r0 = np.dot(r0, r0) ** .5
    sqrt_mu = k**.5
    alpha = -np.dot(v0, v0) / k + 2 / norm_r0

    # First guess
    if alpha > 0:
        # Elliptic orbit
        xi_new = sqrt_mu * tof * alpha
    elif alpha < 0:
        # Hyperbolic orbit
        xi_new = (np.sign(tof) * (-1 / alpha)**.5 *
                  np.log((-2 * k * alpha * tof) /
                         (dot_r0v0 + np.sign(tof) *
                          np.sqrt(-k / alpha) * (1 - norm_r0 * alpha))))
    else:
        # Parabolic orbit
        # (Conservative initial guess)
        xi_new = sqrt_mu * tof / norm_r0

    # Newton-Raphson iteration on the Kepler equation
    count = 0
    while count < numiter:
        xi = xi_new
        psi = xi * xi * alpha
        c2_psi = c2(psi)
        c3_psi = c3(psi)
        norm_r = (xi * xi * c2_psi +
                  dot_r0v0 / sqrt_mu * xi * (1 - psi * c3_psi) +
                  norm_r0 * (1 - psi * c2_psi))
        xi_new = xi + (sqrt_mu * tof - xi * xi * xi * c3_psi -
                       dot_r0v0 / sqrt_mu * xi * xi * c2_psi -
                       norm_r0 * xi * (1 - psi * c3_psi)) / norm_r
        if abs(np.divide(xi_new - xi, xi_new)) < rtol or abs(xi_new - xi) < rtol:
            break
        else:
            count += 1
    else:
        raise RuntimeError("Maximum number of iterations reached")

    # Compute Lagrange coefficients
    f = 1 - xi**2 / norm_r0 * c2_psi
    g = tof - xi**3 / sqrt_mu * c3_psi

    gdot = 1 - xi**2 / norm_r * c2_psi
    fdot = sqrt_mu / (norm_r * norm_r0) * xi * (psi * c3_psi - 1)

    return f, g, fdot, gdot
