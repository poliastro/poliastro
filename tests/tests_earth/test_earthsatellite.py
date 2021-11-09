import numpy as np
import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from poliastro.bodies import Earth, Mars
from poliastro.earth import EarthSatellite
from poliastro.earth.atmosphere import COESA76
from poliastro.earth.enums import EarthGravity
from poliastro.spacecraft import Spacecraft
from poliastro.twobody.orbit import Orbit


def test_earth_satellite_orbit():
    r = [3_539.08827417, 5_310.19903462, 3_066.31301457] * u.km
    v = [-6.49780849, 3.24910291, 1.87521413] * u.km / u.s
    ss = Orbit.from_vectors(Earth, r, v)
    C_D = 2.2 * u.one  # Dimensionless (any value would do)
    A = ((np.pi / 4.0) * (u.m ** 2)).to(u.km ** 2)
    m = 100 * u.kg
    spacecraft = Spacecraft(A, C_D, m)
    earth_satellite = EarthSatellite(ss, spacecraft)
    assert isinstance(earth_satellite.orbit, Orbit)


def test_orbit_attractor():
    r = [3_539.08827417, 5_310.19903462, 3_066.31301457] * u.km
    v = [-6.49780849, 3.24910291, 1.87521413] * u.km / u.s
    ss = Orbit.from_vectors(Mars, r, v)
    C_D = 2.2 * u.one  # Dimensionless (any value would do)
    A = ((np.pi / 4.0) * (u.m ** 2)).to(u.km ** 2)
    m = 100 * u.kg
    spacecraft = Spacecraft(A, C_D, m)
    with pytest.raises(ValueError) as excinfo:
        EarthSatellite(ss, spacecraft)
    assert "The attractor must be Earth" in excinfo.exconly()


def test_propagate_instance():
    tof = 1.0 * u.min
    ss0 = Orbit.from_classical(
        attractor=Earth,
        a=1000 * u.km,
        ecc=0.75 * u.one,
        inc=63.4 * u.deg,
        raan=0 * u.deg,
        argp=270 * u.deg,
        nu=80 * u.deg,
    )
    C_D = 2.2 * u.one  # Dimensionless (any value would do)
    A = ((np.pi / 4.0) * (u.m ** 2)).to(u.km ** 2)
    m = 100 * u.kg
    spacecraft = Spacecraft(A, C_D, m)
    earth_satellite = EarthSatellite(ss0, spacecraft)
    orbit_with_j2 = earth_satellite.propagate(tof=tof, gravity=EarthGravity.J2)
    orbit_without_perturbation = earth_satellite.propagate(tof)
    orbit_with_atmosphere_and_j2 = earth_satellite.propagate(
        tof=tof, gravity=EarthGravity.J2, atmosphere=COESA76()
    )
    assert isinstance(orbit_with_j2, EarthSatellite)
    assert isinstance(orbit_with_atmosphere_and_j2, EarthSatellite)
    assert isinstance(orbit_without_perturbation, EarthSatellite)


# Examples taken from 'Space Mission Analysis and Design', Third Edition, page 156"
@pytest.mark.parametrize(
    "ndays, norbits, inclination, expected_a, expected_ecc",
    [
        (3, 43, 108 * u.deg, 7169 * u.km, 0.004 * u.one),
        (16, 233, 98.2 * u.deg, 7077.7 * u.km, 0.001 * u.one),
        (17, 244, 108.05 * u.deg, 7162.7 * u.km, 0.004 * u.one),
    ],
)
def test_ground_track(ndays, norbits, inclination, expected_a, expected_ecc):
    ss0 = Orbit.from_classical(
        Earth, 5000 * u.km, 0 * u.one, inclination, 0 * u.deg, 0 * u.deg, 0 * u.deg,
    )
    earth_satellite = EarthSatellite(ss0, None)
    orbit = earth_satellite.rgt(ndays, norbits)
    assert_quantity_allclose(orbit.a, expected_a, atol=5e-2 * u.km, rtol=1e-5)
    assert_quantity_allclose(orbit.ecc, expected_ecc, atol=1e-3 * u.one)
