import numpy as np
import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time
from numpy.testing import assert_allclose

from poliastro.bodies import Earth, Moon
from poliastro.maneuver import Maneuver
from poliastro.twobody import Orbit


def test_maneuver_raises_error_if_units_are_wrong():
    wrong_dt = 1.0
    _v = np.zeros(3) * u.km / u.s  # Unused velocity
    with pytest.raises(u.UnitsError) as excinfo:
        Maneuver([wrong_dt, _v])
    assert (
        "UnitsError: Argument 'dts' to function '_initialize' must be in units convertible to 's'."
        in excinfo.exconly()
    )


def test_maneuver_raises_error_if_dvs_are_not_vectors():
    dt = 1 * u.s
    wrong_dv = 1 * u.km / u.s
    with pytest.raises(ValueError) as excinfo:
        Maneuver((dt, wrong_dv))
    assert "ValueError: Delta-V must be three dimensions vectors" in excinfo.exconly()


def test_maneuver_total_time():
    dt1 = 10.0 * u.s
    dt2 = 100.0 * u.s
    _v = np.zeros(3) * u.km / u.s  # Unused velocity
    expected_total_time = 110.0 * u.s
    man = Maneuver((dt1, _v), (dt2, _v))
    assert_quantity_allclose(man.get_total_time(), expected_total_time)


def test_maneuver_impulse():
    dv = [1, 0, 0] * u.m / u.s
    man = Maneuver.impulse(dv)
    assert man.impulses[0] == (0 * u.s, dv)


@pytest.mark.parametrize("nu", [0, 180] * u.deg)
def test_hohmann_maneuver(nu):
    # Data from Vallado, example 6.1
    alt_i = 191.34411 * u.km
    alt_f = 35781.34857 * u.km
    _a = 0 * u.deg
    ss_i = Orbit.from_classical(
        Earth, Earth.R + alt_i, ecc=0 * u.one, inc=_a, raan=_a, argp=_a, nu=nu
    )

    # Expected output
    expected_dv = 3.935224 * u.km / u.s
    expected_t_pericenter = ss_i.time_to_anomaly(0 * u.deg)
    expected_t_trans = 5.256713 * u.h
    expected_total_time = expected_t_pericenter + expected_t_trans

    man = Maneuver.hohmann(ss_i, Earth.R + alt_f)
    assert_quantity_allclose(man.get_total_cost(), expected_dv, rtol=1e-5)
    assert_quantity_allclose(man.get_total_time(), expected_total_time, rtol=1e-5)

    assert_quantity_allclose(
        ss_i.apply_maneuver(man).ecc, 0 * u.one, atol=1e-14 * u.one
    )


@pytest.mark.parametrize("nu", [0, 180] * u.deg)
def test_bielliptic_maneuver(nu):
    # Data from Vallado, example 6.2
    alt_i = 191.34411 * u.km
    alt_b = 503873.0 * u.km
    alt_f = 376310.0 * u.km
    _a = 0 * u.deg
    ss_i = Orbit.from_classical(
        Earth, Earth.R + alt_i, ecc=0 * u.one, inc=_a, raan=_a, argp=_a, nu=nu
    )

    # Expected output
    expected_dv = 3.904057 * u.km / u.s
    expected_t_pericenter = ss_i.time_to_anomaly(0 * u.deg)
    expected_t_trans = 593.919803 * u.h
    expected_total_time = expected_t_pericenter + expected_t_trans

    man = Maneuver.bielliptic(ss_i, Earth.R + alt_b, Earth.R + alt_f)

    assert_allclose(ss_i.apply_maneuver(man).ecc, 0 * u.one, atol=1e-12 * u.one)
    assert_quantity_allclose(man.get_total_cost(), expected_dv, rtol=1e-5)
    assert_quantity_allclose(man.get_total_time(), expected_total_time, rtol=1e-6)


def test_apply_maneuver_correct_dimensions():
    orb = Orbit.from_vectors(
        Moon,
        [-22681.58976181, 942.47776988, 0] * u.km,
        [-0.04578917, -0.19408599, 0.0] * u.km / u.s,
        Time("2023-08-30 23:14", scale="tdb"),
    )
    man = Maneuver((1 * u.s, [0.01, 0, 0] * u.km / u.s))

    new_orb = orb.apply_maneuver(man, intermediate=False)

    assert new_orb.r.ndim == 1
    assert new_orb.v.ndim == 1
