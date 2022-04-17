from functools import partial

import numpy as np
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from hypothesis import example, given, settings, strategies as st

from poliastro.twobody.sampling import sample_closed

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
    result = sample_closed(min_nu, ecc, max_nu)

    assert ((-np.pi * u.rad <= result) & (result <= np.pi * u.rad)).all()


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
def test_sample_closed_starts_at_min_anomaly_if_in_range(min_nu, ecc, max_nu):
    result = sample_closed(min_nu, ecc, max_nu)

    assert_quantity_allclose(result[0], min_nu, atol=1e-15 * u.rad)


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
def test_sample_closed_starts_and_ends_at_min_anomaly_if_in_range_and_no_max_given(
    min_nu, ecc
):
    result = sample_closed(min_nu, ecc)

    assert_quantity_allclose(result[0], min_nu, atol=1e-14 * u.rad)
    assert_quantity_allclose(result[-1], min_nu, atol=1e-14 * u.rad)
