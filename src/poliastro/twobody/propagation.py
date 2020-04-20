"""Propagation algorithms

"""
import functools
from warnings import warn

import numpy as np
from astropy import units as u
from astropy.coordinates import CartesianDifferential, CartesianRepresentation
from scipy.integrate import solve_ivp

from poliastro.core.propagation import (
    func_twobody,
    kepler as kepler_fast,
    mean_motion as mean_motion_fast,
)
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

    rr, vv = method(
        orbit.attractor.k, orbit.r, orbit.v, time_of_flight.to(u.s), rtol=rtol, **kwargs
    )

    # TODO: Turn these into unit tests
    assert rr.ndim == 2
    assert vv.ndim == 2

    cartesian = CartesianRepresentation(
        rr, differentials=CartesianDifferential(vv, xyz_axis=1), xyz_axis=1
    )

    # If the frame supports obstime, set the time values
    kwargs = {}
    if "obstime" in orbit.frame.frame_attributes:
        kwargs["obstime"] = orbit.epoch + time_of_flight
    else:
        warn(
            "Frame {} does not support 'obstime', time values were not returned".format(
                orbit.frame.__class__
            )
        )

    # Use of a protected method instead of frame.realize_frame
    # because the latter does not let the user choose the representation type
    # in one line despite its parameter names, see
    # https://github.com/astropy/astropy/issues/7784
    return orbit.frame._replicate(cartesian, representation_type="cartesian", **kwargs)
