# coding: utf-8
import pytest

import numpy as np
from numpy.testing import assert_almost_equal, assert_array_almost_equal

from astropy import units as u
from astropy.constants import R_earth

from poliastro.bodies import Earth
from poliastro import maneuver


# Object oriented API

def test_maneuver_total_time():
    dt1 = 10.0 * u.s
    dt2 = 100.0 * u.s
    _v = np.zeros(3) * u.km / u.s  # Unused velocity
    expected_total_time = 110.0 * u.s
    man = maneuver.Maneuver([(dt1, _v), (dt2, _v)])
    assert_almost_equal(man.total_time(), expected_total_time)
    assert_array_almost_equal(man.tof, expected_total_time)


# Procedural API

def test_hohmann_example():
    # Data from Vallado, example 6.1
    alt_i = 191.34411 * u.km
    alt_f = 35781.34857 * u.km
    k = Earth.k.to(u.km ** 3 / u.s ** 2)
    R = R_earth.to(u.km)  # TODO: Earth.R.to(u.km)
    expected_dv = 3.935224 * u.km / u.s
    expected_t_trans = 5.256713 * u.h
    r_i = R + alt_i
    r_f = R + alt_f
    dv_a, dv_b, t_trans = maneuver.hohmann(k, r_i, r_f)
    assert_almost_equal(abs(dv_a) + abs(dv_b), expected_dv, decimal=3)
    assert_almost_equal(t_trans.to(u.h), expected_t_trans, decimal=3)


def test_bielliptic_example():
    # Data from Vallado, example 6.2
    alt_i = 191.34411 * u.km
    alt_b = 503873.0 * u.km
    alt_f = 376310.0 * u.km
    k = Earth.k.to(u.km ** 3 / u.s ** 2)
    R = R_earth.to(u.km)
    expected_dv = 3.904057 * u.km / u.s
    expected_t_trans = 593.919803 * u.h
    r_i = R + alt_i
    r_b = R + alt_b
    r_f = R + alt_f
    dva, dvb, dvc, t_trans1, t_trans2 = maneuver.bielliptic(k, r_i, r_b, r_f)
    dv = abs(dva) + abs(dvb) + abs(dvc)
    t_trans = t_trans1 + t_trans2
    assert_almost_equal(dv, expected_dv, decimal=3)
    assert_almost_equal(t_trans.to(u.h), expected_t_trans, decimal=1)
