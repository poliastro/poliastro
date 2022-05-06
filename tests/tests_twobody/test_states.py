from astropy import units as u

from poliastro.bodies import Sun
from poliastro.twobody.states import ClassicalState, RVState


def test_state_has_attractor_given_in_constructor():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    ss = ClassicalState(Sun, _d, _, _a, _a, _a, _a, None)
    assert ss.attractor == Sun


def test_classical_state_has_elements_given_in_constructor():
    # Mars data from HORIZONS at J2000
    a = 1.523679 * u.AU
    ecc = 0.093315 * u.one
    inc = 1.85 * u.deg
    raan = 49.562 * u.deg
    argp = 286.537 * u.deg
    nu = 23.33 * u.deg
    ss = ClassicalState(
        Sun, a * (1 - ecc**2), ecc, inc, raan, argp, nu, None
    )
    assert ss.p == a * (1 - ecc**2)
    assert ss.ecc == ecc
    assert ss.inc == inc
    assert ss.raan == raan
    assert ss.argp == argp
    assert ss.nu == nu


def test_rv_state_has_rv_given_in_constructor():
    r = [1.0, 0.0, 0.0] * u.AU
    v = [0.0, 1.0e-6, 0.0] * u.AU / u.s
    ss = RVState(Sun, r, v, None)
    assert (ss.r == r).all()
    assert (ss.v == v).all()
