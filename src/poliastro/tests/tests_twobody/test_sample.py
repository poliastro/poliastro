import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time

from poliastro.bodies import Earth, Sun
from poliastro.twobody import Orbit
from poliastro.twobody.propagation import cowell, kepler, mean_motion
from poliastro.util import norm


@pytest.mark.parametrize("num_points", [3, 5, 7, 9, 11, 101])
def test_sample_num_points(num_points):
    # Data from Vallado, example 2.4
    r0 = [1131.340, -2282.343, 6672.423] * u.km
    v0 = [-5.64305, 4.30333, 2.42879] * u.km / u.s
    ss0 = Orbit.from_vectors(Earth, r0, v0)

    # TODO: Test against the perigee and apogee
    # expected_ss = ss0.propagate(ss0.period / 2)

    rr = ss0.sample(num_points)

    assert len(rr) == num_points
    # assert_quantity_allclose(rr[num_points // 2].data.xyz, expected_ss.r)


def test_sample_big_orbits():

    # See https://github.com/poliastro/poliastro/issues/265
    ss = Orbit.from_vectors(
        Sun,
        [-9018878.6, -94116055, 22619059] * u.km,
        [-49.950923, -12.948431, -4.2925158] * u.km / u.s,
    )
    positions = ss.sample(15)
    assert len(positions) == 15
