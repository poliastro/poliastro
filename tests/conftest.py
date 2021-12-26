import pytest
from astropy.coordinates import solar_system_ephemeris
from astropy.time import Time

solar_system_ephemeris.set("builtin")


@pytest.fixture
def earth_perihelion():
    return Time("2020-01-05 07:47:00", scale="tdb")
