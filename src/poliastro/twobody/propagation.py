"""Propagation algorithms

"""
import numpy as np
import functools

from scipy.integrate import solve_ivp
from poliastro.integrators import DOP835

from astropy import units as u
from poliastro.twobody.rv import rv2coe
from poliastro.twobody.classical import coe2rv
from poliastro.twobody.angles import nu_to_M, M_to_nu

from poliastro.jit import jit
from poliastro.stumpff import c2, c3


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
    r"""Propagates orbit using mean motion

    .. versionadded:: 0.9.0

    Parameters
    ----------
    orbit : ~poliastro.twobody.orbit.Orbit
        the Orbit object to propagate.
    tof : float
        Time of flight (s).

    Notes
    -----
    This method takes initial :math:`\vec{r}, \vec{v}`, calculates classical orbit parameters,
    increases mean anomaly and performs inverse transformation to get final :math:`\vec{r}, \vec{v}`
    The logic is based on formulae (4), (6) and (7) from http://dx.doi.org/10.1007/s10569-013-9476-9

    """

    k = orbit.attractor.k.to(u.km ** 3 / u.s ** 2).value
    r0 = orbit.r.to(u.km).value
    v0 = orbit.v.to(u.km / u.s).value

    # get the initial true anomaly and orbit parameters that are constant over time
    p, ecc, inc, raan, argp, nu0 = rv2coe(k, r0, v0)

    # get the initial mean anomaly
    M0 = nu_to_M(nu0, ecc)
    # strong elliptic or strong hyperbolic orbits
    if np.abs(ecc - 1.0) > 1e-2:
        a = p / (1.0 - ecc ** 2)
        # given the initial mean anomaly, calculate mean anomaly
        # at the end, mean motion (n) equals sqrt(mu / |a^3|)
        with u.set_enabled_equivalencies(u.dimensionless_angles()):
            M = M0 + tof * np.sqrt(k / np.abs(a ** 3)) * u.rad
            nu = M_to_nu(M, ecc)

    # near-parabolic orbit
    else:
        q = p * np.abs(1.0 - ecc) / np.abs(1.0 - ecc ** 2)
        # mean motion n = sqrt(mu / 2 q^3) for parabolic orbit
        with u.set_enabled_equivalencies(u.dimensionless_angles()):
            M = M0 + tof * np.sqrt(k / 2.0 / (q ** 3))
            nu = M_to_nu(M, ecc)
    with u.set_enabled_equivalencies(u.dimensionless_angles()):
        return coe2rv(k, p, ecc, inc, raan, argp, nu)


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

    f, g, fdot, gdot = _kepler(k, r0, v0, tof, numiter)

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


@jit
def _kepler(k, r0, v0, tof, numiter):
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
        if abs(xi_new - xi) < 1e-7:
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
