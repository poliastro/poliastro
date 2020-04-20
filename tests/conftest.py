import pytest
from astropy.time import Time


@pytest.fixture
def earth_perihelion():
    return Time("2020-01-05 07:47:00", scale="tdb")
