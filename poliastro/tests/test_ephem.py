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


def test_de421_is_selected_if_present():
    kernel_filenames = ["de421.bsp"]
    selected_kernel = ephem.select_kernel(kernel_filenames)
    assert selected_kernel == "de421.bsp"


def test_first_kernel_is_selected_if_de421_not_present():
    kernel_filenames = ["1.bsp", "2.bsp"]
    selected_kernel = ephem.select_kernel(kernel_filenames)
    assert selected_kernel == "1.bsp"


def test_None_is_returned_if_list_is_empty():
    assert ephem.select_kernel([]) is None
