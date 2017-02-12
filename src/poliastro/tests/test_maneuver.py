# coding: utf-8
import pytest

import numpy as np
from numpy.testing import assert_allclose

from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.maneuver import Maneuver


def test_maneuver_raises_error_if_units_are_wrong():
    wrong_dt = 1.0
    _v = np.zeros(3) * u.km / u.s  # Unused velocity
    with pytest.raises(u.UnitsError) as excinfo:
        man = Maneuver([wrong_dt, _v])
    assert ("UnitsError: Argument 'dts' to function '_initialize' must be in units convertible to 's'."
            in excinfo.exconly())


def test_maneuver_raises_error_if_dvs_are_not_vectors():
    dt = 1 * u.s
    wrong_dv = 1 * u.km / u.s
    with pytest.raises(ValueError) as excinfo:
        man = Maneuver((dt, wrong_dv))
    assert ("ValueError: Delta-V must be three dimensions vectors"
            in excinfo.exconly())


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


def test_hohmann_maneuver():
    # Data from Vallado, example 6.1
    alt_i = 191.34411 * u.km
    alt_f = 35781.34857 * u.km
    ss_i = Orbit.circular(Earth, alt_i)
    expected_dv = 3.935224 * u.km / u.s
    expected_t_trans = 5.256713 * u.h

    man = Maneuver.hohmann(ss_i, Earth.R + alt_f)

    assert_quantity_allclose(ss_i.apply_maneuver(man).ecc, 0 * u.one, atol=1e-15 * u.one)
    assert_quantity_allclose(man.get_total_cost(), expected_dv, rtol=1e-6)
    assert_quantity_allclose(man.get_total_time(), expected_t_trans, rtol=1e-6)


def test_bielliptic_maneuver():
    # Data from Vallado, example 6.2
    alt_i = 191.34411 * u.km
    alt_b = 503873.0 * u.km
    alt_f = 376310.0 * u.km
    ss_i = Orbit.circular(Earth, alt_i)
    expected_dv = 3.904057 * u.km / u.s
    expected_t_trans = 593.919803 * u.h

    man = Maneuver.bielliptic(ss_i, Earth.R + alt_b, Earth.R + alt_f)

    assert_allclose(ss_i.apply_maneuver(man).ecc, 0 * u.one, atol=1e-12 * u.one)
    assert_quantity_allclose(man.get_total_cost(), expected_dv, rtol=1e-6)
    assert_quantity_allclose(man.get_total_time(), expected_t_trans, rtol=1e-6)
