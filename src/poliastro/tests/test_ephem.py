# coding: utf-8
from astropy import time
from astropy import units as u

from poliastro import ephem


def test_proper_velocity_units():
    # Bug #49
    _body = "earth"
    _epoch = time.Time("2000-01-01 00:00")

    r, v = ephem.get_body_ephem(_body, _epoch)

    assert r.unit == u.km
    assert v.unit == u.km / u.day
