# coding: utf-8
import pytest

import numpy as np
from numpy.testing import assert_almost_equal

from astropy import units as u

from poliastro.bodies import Earth
from poliastro.twobody import State
from poliastro.maneuver import Maneuver


def test_maneuver_raises_error_if_units_are_wrong():
    wrong_dt = 1.0
    _v = np.zeros(3) * u.km / u.s  # Unused velocity
    with pytest.raises(u.UnitsError) as excinfo:
        man = Maneuver([wrong_dt, _v])
    assert ("UnitsError: Units must be consistent"
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
    assert_almost_equal(man.get_total_time().to(u.s).value,
                        expected_total_time.value)


def test_maneuver_impulse():
    dv = [1, 0, 0] * u.m / u.s
    man = Maneuver.impulse(dv)
    assert man.impulses[0] == (0 * u.s, dv)


def test_hohmann_maneuver():
    # Data from Vallado, example 6.1
    alt_i = 191.34411 * u.km
    alt_f = 35781.34857 * u.km
    ss_i = State.circular(Earth, alt_i)
    expected_dv = 3.935224 * u.km / u.s
    expected_t_trans = 5.256713 * u.h
    man = Maneuver.hohmann(ss_i, Earth.R + alt_f)
    assert_almost_equal(ss_i.apply_maneuver(man).ecc, 0)
    assert_almost_equal(man.get_total_cost().to(u.km / u.s).value,
                        expected_dv.value,
                        decimal=5)
    assert_almost_equal(man.get_total_time().to(u.h).value,
                        expected_t_trans.value,
                        decimal=5)


def test_bielliptic_maneuver():
    # Data from Vallado, example 6.2
    alt_i = 191.34411 * u.km
    alt_b = 503873.0 * u.km
    alt_f = 376310.0 * u.km
    ss_i = State.circular(Earth, alt_i)
    expected_dv = 3.904057 * u.km / u.s
    expected_t_trans = 593.919803 * u.h
    man = Maneuver.bielliptic(ss_i, Earth.R + alt_b, Earth.R + alt_f)
    assert_almost_equal(ss_i.apply_maneuver(man).ecc, 0)
    assert_almost_equal(man.get_total_cost().to(u.km / u.s).value,
                        expected_dv.value,
                        decimal=5)
    assert_almost_equal(man.get_total_time().to(u.h).value,
                        expected_t_trans.value,
                        decimal=2)
