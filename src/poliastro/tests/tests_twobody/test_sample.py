import pytest

from astropy.tests.helper import assert_quantity_allclose
from astropy import units as u
from astropy.time import Time

from poliastro.bodies import Earth
from poliastro.twobody import Orbit


@pytest.mark.parametrize("time_of_flight", [6 * u.h, 2 * u.day])
def test_sample_one_point_equals_propagation_big_deltas(time_of_flight):
    # Data from Vallado, example 2.4
    r0 = [1131.340, -2282.343, 6672.423] * u.km
    v0 = [-5.64305, 4.30333, 2.42879] * u.km / u.s
    ss0 = Orbit.from_vectors(Earth, r0, v0)

    sample_times = Time([ss0.epoch + time_of_flight])

    expected_ss = ss0.propagate(time_of_flight)

    rr = ss0.sample(sample_times)

    assert (rr[0].get_xyz() == expected_ss.r).all()


@pytest.mark.parametrize("time_of_flight", [1 * u.min, 40 * u.min])
def test_sample_one_point_equals_propagation_small_deltas(time_of_flight):
    # FIXME: Time arithmetic is wrong
    # Data from Vallado, example 2.4
    r0 = [1131.340, -2282.343, 6672.423] * u.km
    v0 = [-5.64305, 4.30333, 2.42879] * u.km / u.s
    ss0 = Orbit.from_vectors(Earth, r0, v0)

    sample_times = Time([ss0.epoch + time_of_flight])

    expected_ss = ss0.propagate(time_of_flight)

    rr = ss0.sample(sample_times)

    assert (rr[0].get_xyz() == expected_ss.r).all()
