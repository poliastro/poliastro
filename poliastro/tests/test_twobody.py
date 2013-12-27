import pytest
from astropy import units as u
from astropy import time

from poliastro.bodies import Sun
from poliastro import twobody


def test_state_raises_valueerror_if_meaningless_state():
    _ = 1.0
    with pytest.raises(ValueError) as excinfo:
        ss = twobody.State(Sun, (_,))
    assert ("ValueError: Incorrect number of parameters"
            in excinfo.exconly())


def test_state_has_attractor_given_in_constructor():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    ss = twobody.State(Sun, (_d, _, _a, _a, _a, _a))
    assert ss.attractor == Sun


def test_default_time_for_new_state():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    _body = Sun  # Unused body
    expected_epoch = time.Time("J2000", scale='utc')
    ss = twobody.State(_body, (_d, _, _a, _a, _a, _a))
    assert ss.epoch == expected_epoch


def test_state_has_elements_given_in_constructor():
    # Mars data from HORIZONS at J2000
    a = 1.523679 * u.AU
    ecc = 0.093315
    inc = 1.85 * u.deg
    raan = 49.562 * u.deg
    argp = 286.537 * u.deg
    nu = 23.33 * u.deg
    ss = twobody.State(Sun, (a, ecc, inc, raan, argp, nu))
    assert ss.elements == (a, ecc, inc, raan, argp, nu)


def test_state_raises_unitserror_if_elements_units_are_wrong():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    wrong_angle = 1.0 * u.AU
    with pytest.raises(u.UnitsError) as excinfo:
        ss = twobody.State(Sun, (_d, _, _a, _a, _a, wrong_angle))
    assert ("astropy.units.core.UnitsError: Units must be consistent"
            in excinfo.exconly())


def test_state_raises_valueerror_if_elements_units_are_wrong():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    wrong_angle = 1.0
    with pytest.raises(ValueError) as excinfo:
        ss = twobody.State(Sun, (_d, _, _a, _a, _a, wrong_angle))
    assert ("Elements must have units (use astropy.units)"
            in excinfo.exconly())
