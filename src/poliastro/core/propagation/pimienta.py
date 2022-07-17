import numpy as np
from numba import njit as jit

from poliastro.core.angles import E_to_M, E_to_nu, nu_to_E
from poliastro.core.elements import coe2rv, rv2coe


@jit
def pimienta_coe(k, p, ecc, inc, raan, argp, nu, tof):

    q = p / (1 + ecc)

    # TODO: Do something to allow parabolic and hyperbolic orbits?
    n = np.sqrt(k * (1 - ecc) ** 3 / q**3)
    M0 = E_to_M(nu_to_E(nu, ecc), ecc)

    M = M0 + n * tof

    # Equation (32a), (32b), (32c) and (32d)
    c3 = 5 / 2 + 560 * ecc
    a = 15 * (1 - ecc) / c3
    b = -M / c3
    y = np.sqrt(b**2 / 4 + a**3 / 27)

    # Equation (33)
    x_bar = (-b / 2 + y) ** (1 / 3) - (b / 2 + y) ** (1 / 3)

    # Coefficients from equations (34a) and (34b)
    c15 = 3003 / 14336 + 16384 * ecc
    c13 = 3465 / 13312 - 61440 * ecc
    c11 = 945 / 2816 + 92160 * ecc
    c9 = 175 / 384 - 70400 * ecc
    c7 = 75 / 112 + 28800 * ecc
    c5 = 9 / 8 - 6048 * ecc

    # Precompute x_bar powers, equations (35a) to (35d)
    x_bar2 = x_bar**2
    x_bar3 = x_bar2 * x_bar
    x_bar4 = x_bar3 * x_bar
    x_bar5 = x_bar4 * x_bar
    x_bar6 = x_bar5 * x_bar
    x_bar7 = x_bar6 * x_bar
    x_bar8 = x_bar7 * x_bar
    x_bar9 = x_bar8 * x_bar
    x_bar10 = x_bar9 * x_bar
    x_bar11 = x_bar10 * x_bar
    x_bar12 = x_bar11 * x_bar
    x_bar13 = x_bar12 * x_bar
    x_bar14 = x_bar13 * x_bar
    x_bar15 = x_bar14 * x_bar

    # Function f and its derivatives are given by all the (36) equation set
    f = (
        c15 * x_bar15
        + c13 * x_bar13
        + c11 * x_bar11
        + c9 * x_bar9
        + c7 * x_bar7
        + c5 * x_bar5
        + c3 * x_bar3
        + 15 * (1 - ecc) * x_bar
        - M
    )
    f1 = (
        15 * c15 * x_bar14
        + 13 * c13 * x_bar12
        + 11 * c11 * x_bar10
        + 9 * c9 * x_bar8
        + 7 * c7 * x_bar6
        + 5 * c5 * x_bar4
        + 3 * c3 * x_bar2
        + 15 * (1 - ecc)
    )
    f2 = (
        210 * c15 * x_bar13
        + 156 * c13 * x_bar11
        + 110 * c11 * x_bar9
        + 72 * c9 * x_bar7
        + 42 * c7 * x_bar5
        + 20 * c5 * x_bar3
        + 6 * c3 * x_bar
    )
    f3 = (
        2730 * c15 * x_bar12
        + 1716 * c13 * x_bar10
        + 990 * c11 * x_bar8
        + 504 * c9 * x_bar6
        + 210 * c7 * x_bar4
        + 60 * c5 * x_bar2
        + 6 * c3
    )
    f4 = (
        32760 * c15 * x_bar11
        + 17160 * c13 * x_bar9
        + 7920 * c11 * x_bar7
        + 3024 * c9 * x_bar5
        + 840 * c7 * x_bar3
        + 120 * c5 * x_bar
    )
    f5 = (
        360360 * c15 * x_bar10
        + 154440 * c13 * x_bar8
        + 55440 * c11 * x_bar6
        + 15120 * c9 * x_bar4
        + 2520 * c7 * x_bar2
        + 120 * c5
    )
    f6 = (
        3603600 * c15 * x_bar9
        + 1235520 * c13 * x_bar7
        + 332640 * c11 * x_bar5
        + 60480 * c9 * x_bar3
        + 5040 * c7 * x_bar
    )
    f7 = (
        32432400 * c15 * x_bar8
        + 8648640 * c13 * x_bar6
        + 1663200 * c11 * x_bar4
        + 181440 * c9 * x_bar2
        + 5040 * c7
    )
    f8 = (
        259459200 * c15 * x_bar7
        + 51891840 * c13 * x_bar5
        + 6652800 * c11 * x_bar3
        + 362880 * c9 * x_bar
    )
    f9 = (
        1.8162144e9 * c15 * x_bar6
        + 259459200 * c13 * x_bar4
        + 19958400 * c11 * x_bar2
        + 362880 * c9
    )
    f10 = (
        1.08972864e10 * c15 * x_bar5
        + 1.0378368e9 * c13 * x_bar3
        + 39916800 * c11 * x_bar
    )
    f11 = (
        5.4486432e10 * c15 * x_bar4
        + 3.1135104e9 * c13 * x_bar2
        + 39916800 * c11
    )
    f12 = 2.17945728e11 * c15 * x_bar3 + 6.2270208e9 * c13 * x_bar
    f13 = 6.53837184 * c15 * x_bar2 + 6.2270208e9 * c13
    f14 = 1.307674368e12 * c15 * x_bar
    f15 = 1.307674368e12 * c15

    # Solving g parameters defined by equations (37a), (37b), (37c) and (37d)
    g1 = 1 / 2
    g2 = 1 / 6
    g3 = 1 / 24
    g4 = 1 / 120
    g5 = 1 / 720
    g6 = 1 / 5040
    g7 = 1 / 40320
    g8 = 1 / 362880
    g9 = 1 / 3628800
    g10 = 1 / 39916800
    g11 = 1 / 479001600
    g12 = 1 / 6.2270208e9
    g13 = 1 / 8.71782912e10
    g14 = 1 / 1.307674368e12

    # Solving for the u_{i} and h_{i} variables defined by equation (38)
    u1 = -f / f1

    h2 = f1 + g1 * u1 * f2
    u2 = -f / h2

    h3 = f1 + g1 * u2 * f2 + g2 * u2**2 * f3
    u3 = -f / h3

    h4 = f1 + g1 * u3 * f2 + g2 * u3**2 * f3 + g3 * u3**3 * f4
    u4 = -f / h4

    h5 = (
        f1
        + g1 * u4 * f2
        + g2 * u4**2 * f3
        + g3 * u4**3 * f4
        + g4 * u4**4 * f5
    )
    u5 = -f / h5

    h6 = (
        f1
        + g1 * u5 * f2
        + g2 * u5**2 * f3
        + g3 * u5**3 * f4
        + g4 * u5**4 * f5
        + g5 * u5**5 * f6
    )
    u6 = -f / h6

    h7 = (
        f1
        + g1 * u6 * f2
        + g2 * u6**2 * f3
        + g3 * u6**3 * f4
        + g4 * u6**4 * f5
        + g5 * u6**5 * f6
        + g6 * u6**6 * f7
    )
    u7 = -f / h7

    h8 = (
        f1
        + g1 * u7 * f2
        + g2 * u7**2 * f3
        + g3 * u7**3 * f4
        + g4 * u7**4 * f5
        + g5 * u7**5 * f6
        + g6 * u7**6 * f7
        + g7 * u7**7 * f8
    )
    u8 = -f / h8

    h9 = (
        f1
        + g1 * u8 * f2
        + g2 * u8**2 * f3
        + g3 * u8**3 * f4
        + g4 * u8**4 * f5
        + g5 * u8**5 * f6
        + g6 * u8**6 * f7
        + g7 * u8**7 * f8
        + g8 * u8**8 * f9
    )
    u9 = -f / h9

    h10 = (
        f1
        + g1 * u9 * f2
        + g2 * u9**2 * f3
        + g3 * u9**3 * f4
        + g4 * u9**4 * f5
        + g5 * u9**5 * f6
        + g6 * u9**6 * f7
        + g7 * u9**7 * f8
        + g8 * u9**8 * f9
        + g9 * u9**9 * f10
    )
    u10 = -f / h10

    h11 = (
        f1
        + g1 * u10 * f2
        + g2 * u10**2 * f3
        + g3 * u10**3 * f4
        + g4 * u10**4 * f5
        + g5 * u10**5 * f6
        + g6 * u10**6 * f7
        + g7 * u10**7 * f8
        + g8 * u10**8 * f9
        + g9 * u10**9 * f10
        + g10 * u10**10 * f11
    )
    u11 = -f / h11

    h12 = (
        f1
        + g1 * u11 * f2
        + g2 * u11**2 * f3
        + g3 * u11**3 * f4
        + g4 * u11**4 * f5
        + g5 * u11**5 * f6
        + g6 * u11**6 * f7
        + g7 * u11**7 * f8
        + g8 * u11**8 * f9
        + g9 * u11**9 * f10
        + g10 * u11**10 * f11
        + g11 * u11**11 * f12
    )
    u12 = -f / h12

    h13 = (
        f1
        + g1 * u12 * f2
        + g2 * u12**2 * f3
        + g3 * u12**3 * f4
        + g4 * u12**4 * f5
        + g5 * u12**5 * f6
        + g6 * u12**6 * f7
        + g7 * u12**7 * f8
        + g8 * u12**8 * f9
        + g9 * u12**9 * f10
        + g10 * u12**10 * f11
        + g11 * u12**11 * f12
        + g12 * u12**12 * f13
    )
    u13 = -f / h13

    h14 = (
        f1
        + g1 * u13 * f2
        + g2 * u13**2 * f3
        + g3 * u13**3 * f4
        + g4 * u13**4 * f5
        + g5 * u13**5 * f6
        + g6 * u13**6 * f7
        + g7 * u13**7 * f8
        + g8 * u13**8 * f9
        + g9 * u13**9 * f10
        + g10 * u13**10 * f11
        + g11 * u13**11 * f12
        + g12 * u13**12 * f13
        + g13 * u13**13 * f14
    )
    u14 = -f / h14

    h15 = (
        f1
        + g1 * u14 * f2
        + g2 * u14**2 * f3
        + g3 * u14**3 * f4
        + g4 * u14**4 * f5
        + g5 * u14**5 * f6
        + g6 * u14**6 * f7
        + g7 * u14**7 * f8
        + g8 * u14**8 * f9
        + g9 * u14**9 * f10
        + g10 * u14**10 * f11
        + g11 * u14**11 * f12
        + g12 * u14**12 * f13
        + g13 * u14**13 * f14
        + g14 * u14**14 * f15
    )
    u15 = -f / h15

    # Solving for x
    x = x_bar + u15
    w = x - 0.01171875 * x**17 / (1 + ecc)

    # Solving for the true anomaly from eccentricity anomaly
    E = M + ecc * (
        -16384 * w**15
        + 61440 * w**13
        - 92160 * w**11
        + 70400 * w**9
        - 28800 * w**7
        + 6048 * w**5
        - 560 * w**3
        + 15 * w
    )

    return E_to_nu(E, ecc)


@jit
def pimienta(k, r0, v0, tof):
    """Raw algorithm for Adonis' Pimienta and John L. Crassidis 15th order
    polynomial Kepler solver.

    Parameters
    ----------
    k : float
        Standar Gravitational parameter.
    r0 : numpy.ndarray
        Initial position vector wrt attractor center.
    v0 : numpy.ndarray
        Initial velocity vector.
    tof : float
        Time of flight.

    Returns
    -------
    rr: numpy.ndarray
        Final position vector.
    vv: numpy.ndarray
        Final velocity vector.

    Notes
    -----
    This algorithm was derived from the original paper:
    Pimienta-Peñalver, A. & Crassidis, John. (2013). Accurate Kepler equation
    solver without transcendental function evaluations. Advances in the Astronautical Sciences. 147. 233-247.
    """

    # TODO: implement hyperbolic case

    # Solve first for eccentricity and mean anomaly
    p, ecc, inc, raan, argp, nu = rv2coe(k, r0, v0)
    nu = pimienta_coe(k, p, ecc, inc, raan, argp, nu, tof)

    return coe2rv(k, p, ecc, inc, raan, argp, nu)
