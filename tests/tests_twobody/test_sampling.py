from functools import partial

from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from hypothesis import example, given, settings, strategies as st
import numpy as np
import pytest

from poliastro.examples import iss
from poliastro.twobody.sampling import TrueAnomalyBounds, sample_closed

angles = partial(st.floats, min_value=-2 * np.pi, max_value=2 * np.pi)
eccentricities = partial(st.floats, min_value=0, max_value=1, exclude_max=True)


@st.composite
def with_units(draw, elements, unit):
    angle = draw(elements)
    return angle * unit


angles_q = partial(with_units, elements=angles(), unit=u.rad)
eccentricities_q = partial(with_units, elements=eccentricities(), unit=u.one)


@settings(deadline=None)
@given(
    min_nu=angles_q(),
    ecc=eccentricities_q(),
    max_nu=st.one_of(angles_q(), st.none()),
)
def test_sample_closed_is_always_between_minus_pi_and_pi(min_nu, ecc, max_nu):
    result = sample_closed(ecc, min_nu, max_nu)

    assert ((-np.pi * u.rad <= result) & (result <= np.pi * u.rad)).all()


@pytest.mark.xfail(
    reason="Some corner cases around the circle boundary need closer inspection"
)
@settings(deadline=None)
@given(
    min_nu=with_units(
        elements=st.floats(
            min_value=-np.pi, max_value=np.pi, exclude_max=True
        ),
        unit=u.rad,
    ),
    ecc=eccentricities_q(),
    max_nu=st.one_of(angles_q(), st.none()),
)
@example(0 * u.rad, 0 * u.one, 0 * u.rad)
@example(3.14159265 << u.rad, 0 << u.one, None)
def test_sample_closed_starts_at_min_anomaly_if_in_range(min_nu, ecc, max_nu):
    result = sample_closed(ecc, min_nu, max_nu)

    assert_quantity_allclose(result[0], min_nu, atol=1e-15 * u.rad)


@pytest.mark.xfail(
    reason="Some corner cases around the circle boundary need closer inspection"
)
@settings(deadline=None)
@given(
    min_nu=with_units(
        elements=st.floats(min_value=-np.pi, max_value=np.pi), unit=u.rad
    ),
    ecc=eccentricities_q(),
)
@example(1e-16 * u.rad, 0 * u.one)
@example(0 * u.rad, 0 * u.one)
@example(0 * u.rad, 0.88680956 * u.one)
@example(0 << u.rad, (1 - 1e-16) << u.one)
@example(3.14159265 << u.rad, 0 << u.one)
def test_sample_closed_starts_and_ends_at_min_anomaly_if_in_range_and_no_max_given(
    min_nu, ecc
):
    result = sample_closed(ecc, min_nu)

    # FIXME: In some corner cases the resulting anomaly goes out of range
    # and rather than trying to fix it right now
    # we will inspect this more closely at some point
    assert_quantity_allclose(result[0], min_nu, atol=1e-14 * u.rad)
    assert_quantity_allclose(result[-1], min_nu, atol=1e-7 * u.rad)


@pytest.mark.parametrize("num_values", [3, 5, 7, 9, 11, 101])
def test_sample_num_points(num_values, elliptic):

    # TODO: Test against the perigee and apogee
    # expected_ss = ss0.propagate(ss0.period / 2)

    strategy = TrueAnomalyBounds(num_values=num_values)
    coords, epochs = strategy.sample(elliptic)

    assert len(coords) == len(epochs) == num_values
    # assert_quantity_allclose(rr[num_points // 2].data.xyz, expected_ss.r)


@pytest.mark.parametrize("min_anomaly", [-30 * u.deg, -10 * u.deg])
@pytest.mark.parametrize("max_anomaly", [10 * u.deg, 30 * u.deg])
def test_sample_hyperbolic_limits(hyperbolic, min_anomaly, max_anomaly):
    num_values = 50

    strategy = TrueAnomalyBounds(
        min_nu=min_anomaly, max_nu=max_anomaly, num_values=num_values
    )
    coords, epochs = strategy.sample(hyperbolic)

    assert len(coords) == len(epochs) == num_values


def test_sample_returns_monotonic_increasing_epochs():
    strategy = TrueAnomalyBounds(num_values=10)
    coords, epochs = strategy.sample(iss)

    assert (np.diff(epochs.jd) > 0).all()
