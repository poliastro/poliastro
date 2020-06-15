import numpy as np
import pytest
from astropy import units as u

from poliastro.atmosphere import COESA76
from poliastro.bodies import Earth, Mars
from poliastro.earth.__init__ import EarthSatellite
from poliastro.earth.enums import EarthGravity
from poliastro.twobody.orbit import Orbit


def test_earth_satellite_orbit():
    r = [3_539.08827417, 5_310.19903462, 3_066.31301457] * u.km
    v = [-6.49780849, 3.24910291, 1.87521413] * u.km / u.s
    ss = Orbit.from_vectors(Earth, r, v)
    earth_satellite = EarthSatellite(ss)
    assert isinstance(earth_satellite.orbit, Orbit)


def test_orbit_attractor():
    r = [3_539.08827417, 5_310.19903462, 3_066.31301457] * u.km
    v = [-6.49780849, 3.24910291, 1.87521413] * u.km / u.s
    ss = Orbit.from_vectors(Mars, r, v)
    with pytest.raises(ValueError) as excinfo:
        EarthSatellite(ss)
    assert "The attractor must be Earth" in excinfo.exconly()


def test_from_vectors():
    r = [1_131.340, -2_282.343, 6_672.423] * u.km
    v = [-5.64305, 4.30333, 2.42879] * u.km / u.s
    earth_sat = EarthSatellite.from_vectors(r, v)
    assert isinstance(earth_sat, EarthSatellite)
    assert (earth_sat.orbit.r == r).all()
    assert (earth_sat.orbit.v == v).all()


def test_from_classical_instance():
    earth_sat = EarthSatellite.from_classical(
        26600 * u.km, 0.75 * u.one, 63.4 * u.deg, 0 * u.deg, 270 * u.deg, 80 * u.deg
    )
    assert isinstance(earth_sat, EarthSatellite)
    assert earth_sat.orbit.a == 26600 * u.km


def test_propagate_instance():
    tof = 1.0 * u.min
    ss0 = EarthSatellite.from_classical(
        1000 * u.km, 0.75 * u.one, 63.4 * u.deg, 0 * u.deg, 270 * u.deg, 80 * u.deg
    )
    orbit_with_j2 = ss0.propagate(tof=tof, gravity=EarthGravity.J2)
    A_over_m = ((np.pi / 4.0) * (u.m ** 2) / (100 * u.kg)).to_value(u.km ** 2 / u.kg)
    orbit_without_perturbation = ss0.propagate(tof)
    orbit_with_atmosphere_and_j2 = ss0.propagate(
        tof=tof, gravity=EarthGravity.J2, atmosphere=COESA76(), A_over_m=A_over_m
    )
    assert isinstance(orbit_with_j2, EarthSatellite)
    assert isinstance(orbit_with_atmosphere_and_j2, EarthSatellite)
    assert isinstance(orbit_without_perturbation, EarthSatellite)
