from astropy import units as u
from astropy.coordinates import solar_system_ephemeris
from astropy.time import Time
import pytest

from poliastro.bodies import Earth, Sun
from poliastro.twobody import Orbit

solar_system_ephemeris.set("builtin")


@pytest.fixture
def earth_perihelion():
    return Time("2020-01-05 07:47:00", scale="tdb")


@pytest.fixture
def elliptic():
    # Data from Vallado, example 2.4
    r0 = [1_131.340, -2_282.343, 6_672.423] << u.km
    v0 = [-5.64305, 4.30333, 2.42879] << (u.km / u.s)
    return Orbit.from_vectors(Earth, r0, v0)


@pytest.fixture()
def hyperbolic():
    r = [
        1.197659243752796e09,
        -4.443716685978071e09,
        -1.747610548576734e09,
    ] * u.km
    v = (
        [5.540549267188614e00, -1.251544669134140e01, -4.848892572767733e00]
        * u.km
        / u.s
    )
    epoch = Time("2015-07-14 07:59", scale="tdb")
    return Orbit.from_vectors(Sun, r, v, epoch)
