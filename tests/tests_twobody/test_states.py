import numpy as np
import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from poliastro.bodies import Earth, Sun
from poliastro.twobody.states import ClassicalState, RVState


def test_state_has_attractor_given_in_constructor():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    ss = ClassicalState(Sun, (_d, _, _a, _a, _a, _a), None)
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
        Sun, (a * (1 - ecc**2), ecc, inc, raan, argp, nu), None
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
    ss = RVState(Sun, (r, v), None)
    assert (ss.r == r).all()
    assert (ss.v == v).all()


def test_mean_motion():
    # From Vallado Example 1-1.
    attractor = Earth
    period = 86164.090518 * u.s
    a = 42164.1696 * u.km
    # Unused variables.
    _ecc = 0 * u.one
    _inc = 1.85 * u.deg
    _raan = 50 * u.deg
    _argp = 200 * u.deg
    _nu = 20 * u.deg

    ss = ClassicalState(
        attractor, (a * (1 - _ecc**2), _ecc, _inc, _raan, _argp, _nu), None
    )

    expected_mean_motion = (2 * np.pi / period) * u.rad
    n = ss.n

    assert_quantity_allclose(n, expected_mean_motion)


def test_coe_to_mee_raises_singularity_error_orbit_equatorial_and_retrograde():
    a = 10000 * u.km
    ecc = 0.3 * u.one
    inc = 180 * u.deg  # True retrograde equatorial case.
    raan = 49.562 * u.deg
    argp = 286.537 * u.deg
    nu = 23.33 * u.deg

    ss = ClassicalState(
        Sun, (a * (1 - ecc**2), ecc, inc, raan, argp, nu), None
    )
    with pytest.raises(ValueError) as excinfo:
        ss.to_equinoctial()
    assert (
        "Cannot compute modified equinoctial set for 180 degrees orbit inclination due to `h` and `k` singularity."
        in excinfo.exconly()
    )
