# coding: utf-8
import numpy as np

from astropy import time
from astropy import units as u

from poliastro import ephem


class FakeKernel(object):
    def __getitem__(self, index):
        return FakeSegment()


class FakeSegment(object):
    def compute_and_differentiate(self, jd1, jd2=None):
        r = np.array([1, 1, 1])
        v = np.array([1, 1, 1])
        return r, v


def test_proper_velocity_units():
    # Bug #49
    _body = 0
    _epoch = time.Time("2000-01-01 00:00")

    r, v = ephem.planet_ephem(_body, _epoch, FakeKernel())

    assert v.unit == u.km / u.day
