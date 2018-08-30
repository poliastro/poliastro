import pytest

from astropy.tests.helper import assert_quantity_allclose
from astropy import units as u
from astropy.time import Time

from poliastro.bodies import Earth, Sun
from poliastro.twobody import Orbit
from poliastro.twobody.propagation import kepler, mean_motion, cowell

from poliastro.util import norm


def test_sample_angle_zero_returns_same():
    # Data from Vallado, example 2.4
    r0 = [1131.340, -2282.343, 6672.423] * u.km
    v0 = [-5.64305, 4.30333, 2.42879] * u.km / u.s
    ss0 = Orbit.from_vectors(Earth, r0, v0)

    nu_values = [0] * u.deg
    _, rr = ss0.sample(ss0.nu + nu_values)

    assert_quantity_allclose(rr[0].get_xyz(), ss0.r)


@pytest.mark.parametrize("time_of_flight", [1 * u.min, 40 * u.min])
@pytest.mark.parametrize("method", [kepler, mean_motion, cowell])
def test_sample_one_point_equals_propagation_small_deltas(time_of_flight, method):
    # Time arithmetic loses precision, see
    # https://github.com/astropy/astropy/issues/6638
    # Data from Vallado, example 2.4
    r0 = [1131.340, -2282.343, 6672.423] * u.km
    v0 = [-5.64305, 4.30333, 2.42879] * u.km / u.s
    ss0 = Orbit.from_vectors(Earth, r0, v0)

    sample_times = Time([ss0.epoch + time_of_flight])

    expected_ss = ss0.propagate(time_of_flight, method)

    _, rr = ss0.sample(sample_times, method)

    assert_quantity_allclose(rr[0].get_xyz(), expected_ss.r)


def test_sampling_cowell_dense_output_agrees_mean_motion():
    r0 = [1131.340, -2282.343, 6672.423] * u.km
    v0 = [-5.64305, 4.30333, 2.42879] * u.km / u.s
    ss0 = Orbit.from_vectors(Earth, r0, v0)
    sample_times = Time([ss0.epoch + x * u.h for x in range(100)])

    _, rr = ss0.sample(sample_times, cowell)
    for i, sample_time in enumerate(sample_times):
        expected_ss = ss0.propagate(sample_time - ss0.epoch)
        assert_quantity_allclose(rr[i].get_xyz(), expected_ss.r)


@pytest.mark.parametrize("time_of_flight", [6 * u.h, 2 * u.day])
@pytest.mark.parametrize("method", [
    kepler,
    mean_motion,
    cowell])
def test_sample_one_point_equals_propagation_big_deltas(time_of_flight, method):
    # Data from Vallado, example 2.4
    r0 = [1131.340, -2282.343, 6672.423] * u.km
    v0 = [-5.64305, 4.30333, 2.42879] * u.km / u.s
    ss0 = Orbit.from_vectors(Earth, r0, v0)

    sample_times = Time([ss0.epoch + time_of_flight])

    expected_ss = ss0.propagate(time_of_flight)

    _, rr = ss0.sample(sample_times, method)

    assert_quantity_allclose(rr[0].get_xyz(), expected_ss.r)


def test_sample_nu_values():
    # Data from Vallado, example 2.4
    r0 = [1131.340, -2282.343, 6672.423] * u.km
    v0 = [-5.64305, 4.30333, 2.42879] * u.km / u.s
    ss0 = Orbit.from_vectors(Earth, r0, v0)

    nu_values = [0, 90, 180] * u.deg

    expected_ss = ss0.propagate(ss0.period / 2)

    _, rr = ss0.sample(nu_values)

    assert len(rr) == len(nu_values)
    assert_quantity_allclose(norm(rr[0].get_xyz()), expected_ss.r_p)
    assert_quantity_allclose(norm(rr[-1].get_xyz()), expected_ss.r_a)


@pytest.mark.parametrize("num_points", [3, 5, 7, 9, 11, 101])
def test_sample_num_points(num_points):
    # Data from Vallado, example 2.4
    r0 = [1131.340, -2282.343, 6672.423] * u.km
    v0 = [-5.64305, 4.30333, 2.42879] * u.km / u.s
    ss0 = Orbit.from_vectors(Earth, r0, v0)

    # TODO: Test against the perigee and apogee
    # expected_ss = ss0.propagate(ss0.period / 2)

    _, rr = ss0.sample(num_points)

    assert len(rr) == num_points
    # assert_quantity_allclose(rr[num_points // 2].data.xyz, expected_ss.r)


@pytest.mark.parametrize('method', [
    mean_motion,
    cowell,
    kepler,
])
def test_sample_big_orbits(method):
    # See https://github.com/poliastro/poliastro/issues/265
    ss = Orbit.from_vectors(
        Sun,
        [-9018878.6, -94116055, 22619059] * u.km,
        [-49.950923, -12.948431, -4.2925158] * u.km / u.s
    )
    times, positions = ss.sample(15, method=method)
    assert len(times) == len(positions) == 15
