"""The following script holds the different high level functions for the
different propagators available at poliastro:

+-------------+------------+-----------------+-----------------+
|  Propagator | Elliptical |    Parabolic    |    Hyperbolic   |
+-------------+------------+-----------------+-----------------+
|  farnocchia |      ✓     |        ✓        |        ✓        |
+-------------+------------+-----------------+-----------------+
|   vallado   |      ✓     |        ✓        |        ✓        |
+-------------+------------+-----------------+-----------------+
|   mikkola   |      ✓     | NOT IMPLEMENTED |        ✓        |
+-------------+------------+-----------------+-----------------+
|   markley   |      ✓     |        x        |        x        |
+-------------+------------+-----------------+-----------------+
|   pimienta  |      ✓     |        ✓        | NOT IMPLEMENTED |
+-------------+------------+-----------------+-----------------+
|   gooding   |      ✓     | NOT IMPLEMENTED | NOT IMPLEMENTED |
+-------------+------------+-----------------+-----------------+
|    danby    |      ✓     |        x        |        ✓        |
+-------------+------------+-----------------+-----------------+
|    cowell   |      ✓     |        ✓        |        ✓        |
+-------------+------------+-----------------+-----------------+

"""
import numpy as np
from astropy import units as u
from astropy.coordinates import CartesianDifferential, CartesianRepresentation
from scipy.integrate import DOP853, solve_ivp

from poliastro.core.propagation import (
    danby as danby_fast,
    farnocchia as farnocchia_fast,
    func_twobody,
    gooding as gooding_fast,
    markley as markley_fast,
    mikkola as mikkola_fast,
    pimienta as pimienta_fast,
    vallado as vallado_fast,
)


def cowell(k, r, v, tofs, rtol=1e-11, *, events=None, f=func_twobody):
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
    events : function(t, u(t)), optional
        passed to solve_ivp: integration stops when this function
        returns <= 0., assuming you set events.terminal=True
    f : function(t0, u, k), optional
        Objective function, default to Keplerian-only forces.

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

    result = solve_ivp(
        f,
        (0, max(tofs)),
        u0,
        args=(k,),
        rtol=rtol,
        atol=1e-12,
        method=DOP853,
        dense_output=True,
        events=events,
    )
    if not result.success:
        raise RuntimeError("Integration failed")

    t_end = (
        min(result.t_events[0]) if result.t_events and len(result.t_events[0]) else None
    )

    rrs = []
    vvs = []
    for i in range(len(tofs)):
        t = tofs[i]
        if t_end is not None and t > t_end:
            t = t_end
        y = result.sol(t)
        rrs.append(y[:3])
        vvs.append(y[3:])

    return rrs * u.km, vvs * u.km / u.s


def farnocchia(k, r, v, tofs, **kwargs):
    """Propagates orbit.

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

    results = [farnocchia_fast(k, r0, v0, tof) for tof in tofs]
    # TODO: Rewrite to avoid iterating twice
    return (
        [result[0] for result in results] * u.km,
        [result[1] for result in results] * u.km / u.s,
    )


def vallado(k, r, v, tofs, numiter=350, **kwargs):
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
    f, g, fdot, gdot = vallado_fast(k, r0, v0, tof, numiter)

    assert np.abs(f * gdot - fdot * g - 1) < 1e-5  # Fixed tolerance

    # Return position and velocity vectors
    r = f * r0 + g * v0
    v = fdot * r0 + gdot * v0

    return r, v


def mikkola(k, r, v, tofs, rtol=None):
    """Solves Kepler Equation by a cubic approximation. This method is valid
    no mater the orbit's nature.

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
    rtol: float
        This method does not require of tolerance since it is non iterative.

    Returns
    -------
    rr : ~astropy.units.Quantity
        Propagated position vectors.
    vv : ~astropy.units.Quantity

    Note
    ----
    This method was derived by Seppo Mikola in his paper *A Cubic Approximation
    For Kepler's Equation* with DOI: https://doi.org/10.1007/BF01235850
    """

    k = k.to(u.m ** 3 / u.s ** 2).value
    r0 = r.to(u.m).value
    v0 = v.to(u.m / u.s).value
    tofs = tofs.to(u.s).value

    results = [mikkola_fast(k, r0, v0, tof) for tof in tofs]
    return (
        [result[0] for result in results] * u.m,
        [result[1] for result in results] * u.m / u.s,
    )


def markley(k, r, v, tofs, rtol=None):
    """Elliptical Kepler Equation solver based on a fifth-order
    refinement of the solution of a cubic equation.

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
    rtol: float
        This method does not require of tolerance since it is non iterative.

    Returns
    -------
    rr : ~astropy.units.Quantity
        Propagated position vectors.
    vv : ~astropy.units.Quantity
        Propagated velocity vectors.

    Note
    ----
    This method was originally presented by Markley in his paper *Kepler Equation Solver*
    with DOI: https://doi.org/10.1007/BF00691917
    """

    k = k.to(u.m ** 3 / u.s ** 2).value
    r0 = r.to(u.m).value
    v0 = v.to(u.m / u.s).value
    tofs = tofs.to(u.s).value

    results = [markley_fast(k, r0, v0, tof) for tof in tofs]
    return (
        [result[0] for result in results] * u.m,
        [result[1] for result in results] * u.m / u.s,
    )


def pimienta(k, r, v, tofs, rtol=None):
    """Kepler solver for both elliptic and parabolic orbits based on a 15th
    order polynomial with accuracies around 10e-5 for elliptic case and 10e-13
    in the hyperbolic regime.

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
    rtol: float
        This method does not require of tolerance since it is non iterative.

    Returns
    -------
    rr : ~astropy.units.Quantity
        Propagated position vectors.
    vv : ~astropy.units.Quantity
        Propagated velocity vectors.

    Note
    ----
    This algorithm was developed by Pimienta-Peñalver and John L. Crassidis in
    their paper *Accurate Kepler Equation solver without trascendental function
    evaluations*. Original paper is on Buffalo's UBIR repository: http://hdl.handle.net/10477/50522
    """

    k = k.to(u.m ** 3 / u.s ** 2).value
    r0 = r.to(u.m).value
    v0 = v.to(u.m / u.s).value
    tofs = tofs.to(u.s).value

    results = [pimienta_fast(k, r0, v0, tof) for tof in tofs]
    return (
        [result[0] for result in results] * u.m,
        [result[1] for result in results] * u.m / u.s,
    )


def gooding(k, r, v, tofs, numiter=150, rtol=1e-8):
    """Solves the Elliptic Kepler Equation with a cubic convergence and
    accuracy better than 10e-12 rad is normally achieved. It is not valid for
    eccentricities equal or greater than 1.0.

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
    rtol: float
        This method does not require of tolerance since it is non iterative.

    Returns
    -------
    rr : ~astropy.units.Quantity
        Propagated position vectors.
    vv : ~astropy.units.Quantity

    Note
    ----
    This method was developed by Gooding and Odell in their paper *The
    hyperbolic Kepler equation (and the elliptic equation revisited)* with
    DOI: https://doi.org/10.1007/BF01235540
    """

    k = k.to(u.m ** 3 / u.s ** 2).value
    r0 = r.to(u.m).value
    v0 = v.to(u.m / u.s).value
    tofs = tofs.to(u.s).value

    results = [gooding_fast(k, r0, v0, tof, numiter=numiter, rtol=rtol) for tof in tofs]
    return (
        [result[0] for result in results] * u.m,
        [result[1] for result in results] * u.m / u.s,
    )


def danby(k, r, v, tofs, rtol=1e-8):
    """Kepler solver for both elliptic and parabolic orbits based on Danby's
    algorithm.

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
    rtol: float
        Relative error for accuracy of the method.

    Returns
    -------
    rr : ~astropy.units.Quantity
        Propagated position vectors.
    vv : ~astropy.units.Quantity
        Propagated velocity vectors.

    Note
    ----
    This algorithm was developed by Danby in his paper *The solution of Kepler
    Equation* with DOI: https://doi.org/10.1007/BF01686811
    """

    k = k.to(u.m ** 3 / u.s ** 2).value
    r0 = r.to(u.m).value
    v0 = v.to(u.m / u.s).value
    tofs = tofs.to(u.s).value

    results = [danby_fast(k, r0, v0, tof) for tof in tofs]
    return (
        [result[0] for result in results] * u.m,
        [result[1] for result in results] * u.m / u.s,
    )


def propagate(orbit, time_of_flight, *, method=farnocchia, rtol=1e-10, **kwargs):
    """Propagate an orbit some time and return the result.

    Parameters
    ----------
    orbit : ~poliastro.twobody.Orbit
        Orbit object to propagate.
    time_of_flight : ~astropy.time.TimeDelta
        Time of propagation.
    method : callable, optional
        Propagation method, default to farnocchia.
    rtol : float, optional
        Relative tolerance, default to 1e-10.

    Returns
    -------
    astropy.coordinates.CartesianRepresentation
        Propagation coordinates.

    """

    # Check if propagator fulfills orbit requirements
    if orbit.ecc < 1.0 and method not in ELLIPTIC_PROPAGATORS:
        raise ValueError(
            "Can not use an parabolic/hyperbolic propagator for elliptical orbits."
        )
    elif orbit.ecc == 1.0 and method not in PARABOLIC_PROPAGATORS:
        raise ValueError(
            "Can not use an elliptic/hyperbolic propagator for parabolic orbits."
        )
    elif orbit.ecc > 1.0 and method not in HYPERBOLIC_PROPAGATORS:
        raise ValueError(
            "Can not use an elliptic/parabolic propagator for hyperbolic orbits."
        )
    else:
        pass

    rr, vv = method(
        orbit.attractor.k,
        orbit.r,
        orbit.v,
        time_of_flight.reshape(-1).to(u.s),
        rtol=rtol,
        **kwargs
    )

    # TODO: Turn these into unit tests
    assert rr.ndim == 2
    assert vv.ndim == 2

    cartesian = CartesianRepresentation(
        rr, differentials=CartesianDifferential(vv, xyz_axis=1), xyz_axis=1
    )

    return cartesian


ELLIPTIC_PROPAGATORS = [
    farnocchia,
    vallado,
    mikkola,
    markley,
    pimienta,
    gooding,
    danby,
    cowell,
]
PARABOLIC_PROPAGATORS = [farnocchia, vallado, mikkola, pimienta, gooding, cowell]
HYPERBOLIC_PROPAGATORS = [
    farnocchia,
    vallado,
    mikkola,
    pimienta,
    gooding,
    danby,
    cowell,
]
ALL_PROPAGATORS = list(
    set(ELLIPTIC_PROPAGATORS) & set(PARABOLIC_PROPAGATORS) & set(HYPERBOLIC_PROPAGATORS)
)
