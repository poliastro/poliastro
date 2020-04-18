import numpy as np

from .._jit import jit
from ..elements import coe2rv, rv2coe
from ..angles import nu_to_M, M_to_nu


@jit
def farnocchia(k, r0, v0, tof):
    r"""Propagates orbit using mean motion.

    This algorithm depends on the geometric shape of the orbit.
    For the case of the strong elliptic or strong hyperbolic orbits:

    ..  math::

        M = M_{0} + \frac{\mu^{2}}{h^{3}}\left ( 1 -e^{2}\right )^{\frac{3}{2}}t

    .. versionadded:: 0.9.0

    Parameters
    ----------
    k : float
        Standar Gravitational parameter
    r0 : ~astropy.units.Quantity
        Initial position vector wrt attractor center.
    v0 : ~astropy.units.Quantity
        Initial velocity vector.
    tof : float
        Time of flight (s).

    Note
    ----
    This method takes initial :math:`\vec{r}, \vec{v}`, calculates classical orbit parameters,
    increases mean anomaly and performs inverse transformation to get final :math:`\vec{r}, \vec{v}`
    The logic is based on formulae (4), (6) and (7) from http://dx.doi.org/10.1007/s10569-013-9476-9

    """

    # get the initial true anomaly and orbit parameters that are constant over time
    p, ecc, inc, raan, argp, nu0 = rv2coe(k, r0, v0)

    # get the initial mean anomaly
    M0 = nu_to_M(nu0, ecc)
    if np.abs(ecc - 1.0) > 1e-2:
        # strong elliptic or strong hyperbolic orbits
        a = p / (1.0 - ecc ** 2)
        n = np.sqrt(k / np.abs(a ** 3))
    else:
        # near-parabolic orbit
        q = p / np.abs(1.0 + ecc)
        n = np.sqrt(k / 2.0 / (q ** 3))

    M = M0 + tof * n
    nu = M_to_nu(M, ecc)

    return coe2rv(k, p, ecc, inc, raan, argp, nu)
