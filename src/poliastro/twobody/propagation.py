"""Propagation algorithms

"""
import numpy as np
import functools

from scipy.integrate import solve_ivp

from poliastro.core.propagation import (
    mean_motion as mean_motion_fast,
    kepler as kepler_fast
)
from poliastro.integrators import DOP835

from astropy import units as u


def func_twobody(t0, u_, k, ad, ad_kwargs):
    """Differential equation for the initial value two body problem.

    This function follows Cowell's formulation.

    Parameters
    ----------
    t0 : float
        Time.
    u_ : ~numpy.ndarray
        Six component state vector [x, y, z, vx, vy, vz] (km, km/s).
    k : float
        Standard gravitational parameter.
    ad : function(t0, u, k)
        Non Keplerian acceleration (km/s2).
    ad_kwargs : optional
        perturbation parameters passed to ad
    """
    ax, ay, az = ad(t0, u_, k, **ad_kwargs)

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


def cowell(orbit, tof, rtol=1e-11, *, ad=None, **ad_kwargs):
    """Propagates orbit using Cowell's formulation.

    Parameters
    ----------
    orbit : ~poliastro.twobody.orbit.Orbit
        the Orbit object to propagate.
    ad : function(t0, u, k), optional
         Non Keplerian acceleration (km/s2), default to None.
    tof : Multiple options
        Time to propagate, float (s),
        Times to propagate, array of float (s).
    rtol : float, optional
        Maximum relative error permitted, default to 1e-10.

    Raises
    ------
    RuntimeError
        If the algorithm didn't converge.

    Note
    -----
    This method uses a Dormand & Prince method of order 8(5,3) available
    in the :py:class:`poliastro.integrators` module. If multiple tofs
    are provided, the method propagates to the maximum value and
    calculates the orher values via dense output

    """
    k = orbit.attractor.k.to(u.km ** 3 / u.s ** 2).value
    x, y, z = orbit.r.to(u.km).value
    vx, vy, vz = orbit.v.to(u.km / u.s).value

    u0 = np.array([x, y, z, vx, vy, vz])

    # Set the non Keplerian acceleration
    if ad is None:
        ad = lambda t0, u_, k_: (0, 0, 0)

    f_with_ad = functools.partial(func_twobody, k=k, ad=ad, ad_kwargs=ad_kwargs)

    multiple_input = hasattr(tof, "__len__")
    if not multiple_input:
        tof = [tof]

    result = solve_ivp(f_with_ad, (0, max(tof)), u0,
                       rtol=rtol, atol=1e-12, method=DOP835,
                       dense_output=True)
    if not result.success:
        raise RuntimeError("Integration failed")

    rrs = []
    vvs = []
    for i in range(len(tof)):
        y = result.sol(tof[i])
        rrs.append(y[:3])
        vvs.append(y[3:])

    if not multiple_input:
        return rrs[0], vvs[0]
    return rrs, vvs


def mean_motion(orbit, tof, **kwargs):
    k = orbit.attractor.k.to(u.km ** 3 / u.s ** 2).value
    r0 = orbit.r.to(u.km).value
    v0 = orbit.v.to(u.km / u.s).value

    return mean_motion_fast(k, r0, v0, tof)


def kepler(orbit, tof, *, numiter=350, **kwargs):
    """Propagates Keplerian orbit.

    Parameters
    ----------
    orbit : ~poliastro.twobody.orbit.Orbit
        the Orbit object to propagate.
    tof : float
        Time of flight (s).
    numiter : int, optional
        Maximum number of iterations, default to 35.

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
    # Compute Lagrange coefficients
    k = orbit.attractor.k.to(u.km ** 3 / u.s ** 2).value
    r0 = orbit.r.to(u.km).value
    v0 = orbit.v.to(u.km / u.s).value

    f, g, fdot, gdot = kepler_fast(k, r0, v0, tof, numiter)

    assert np.abs(f * gdot - fdot * g - 1) < 1e-5  # Fixed tolerance

    # Return position and velocity vectors
    r = f * r0 + g * v0
    v = fdot * r0 + gdot * v0

    return r, v


def propagate(orbit, time_of_flight, *, method=mean_motion, rtol=1e-10, **kwargs):
    """Propagate an orbit some time and return the result.

    """
    r, v = method(orbit, time_of_flight.to(u.s).value, rtol=rtol, **kwargs)
    return orbit.from_vectors(orbit.attractor, r * u.km, v * u.km / u.s, orbit.epoch + time_of_flight)
