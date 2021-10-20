import warnings

import numpy as np
import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time
from numpy.testing import assert_allclose

from poliastro.bodies import Earth, Mercury, Moon
from poliastro.maneuver import Maneuver
from poliastro.twobody import Orbit


def test_maneuver_constructor_raises_error_if_invalid_delta_v():
    dv1 = np.zeros(3) * u.km / u.s
    dv2 = np.ones(2) * u.km / u.s  # Incorrect dv
    with pytest.raises(ValueError) as excinfo:
        with warnings.catch_warnings():
            # Different length numpy arrays generate a deprecation warning.
            warnings.simplefilter("ignore", category=np.VisibleDeprecationWarning)
            Maneuver((0 * u.s, dv1), (2 * u.s, dv2))
    assert "Delta-V must be three dimensions vectors" in excinfo.exconly()


def test_maneuver_raises_error_if_units_are_wrong():
    wrong_dt = 1.0
    _v = np.zeros(3) * u.km / u.s  # Unused velocity
    with pytest.raises(u.UnitsError) as excinfo:
        Maneuver([wrong_dt, _v])
    assert (
        "Argument 'dts' to function '_initialize' must be in units convertible to 's'."
        in excinfo.exconly()
    )


def test_maneuver_raises_error_if_dvs_are_not_vectors():
    dt = 1 * u.s
    wrong_dv = 1 * u.km / u.s
    with pytest.raises(ValueError) as excinfo:
        Maneuver((dt, wrong_dv))
    assert "Delta-V must be three dimensions vectors" in excinfo.exconly()


def test_maneuver_total_time():
    dt1 = 10.0 * u.s
    dt2 = 100.0 * u.s
    _v = np.zeros(3) * u.km / u.s  # Unused velocity
    expected_total_time = 110.0 * u.s
    man = Maneuver((dt1, _v), (dt2, _v))
    assert_quantity_allclose(man.get_total_time(), expected_total_time)


def test_maneuver_impulse():
    dv = [1, 0, 0] * u.m / u.s
    man = Maneuver.impulse(dv)
    assert man.impulses[0] == (0 * u.s, dv)


@pytest.mark.parametrize("nu", [0, -180] * u.deg)
def test_hohmann_maneuver(nu):
    # Data from Vallado, example 6.1
    alt_i = 191.34411 * u.km
    alt_f = 35781.34857 * u.km
    _a = 0 * u.deg
    ss_i = Orbit.from_classical(
        Earth, Earth.R + alt_i, ecc=0 * u.one, inc=_a, raan=_a, argp=_a, nu=nu
    )

    # Expected output
    expected_dv = 3.935224 * u.km / u.s
    expected_t_pericenter = ss_i.time_to_anomaly(0 * u.deg)
    expected_t_trans = 5.256713 * u.h
    expected_total_time = expected_t_pericenter + expected_t_trans

    man = Maneuver.hohmann(ss_i, Earth.R + alt_f)
    assert_quantity_allclose(man.get_total_cost(), expected_dv, rtol=1e-5)
    assert_quantity_allclose(man.get_total_time(), expected_total_time, rtol=1e-5)

    assert_quantity_allclose(
        ss_i.apply_maneuver(man).ecc, 0 * u.one, atol=1e-14 * u.one
    )


@pytest.mark.parametrize("nu", [0, -180] * u.deg)
def test_bielliptic_maneuver(nu):
    # Data from Vallado, example 6.2
    alt_i = 191.34411 * u.km
    alt_b = 503873.0 * u.km
    alt_f = 376310.0 * u.km
    _a = 0 * u.deg
    ss_i = Orbit.from_classical(
        Earth, Earth.R + alt_i, ecc=0 * u.one, inc=_a, raan=_a, argp=_a, nu=nu
    )

    # Expected output
    expected_dv = 3.904057 * u.km / u.s
    expected_t_pericenter = ss_i.time_to_anomaly(0 * u.deg)
    expected_t_trans = 593.919803 * u.h
    expected_total_time = expected_t_pericenter + expected_t_trans

    man = Maneuver.bielliptic(ss_i, Earth.R + alt_b, Earth.R + alt_f)

    assert_allclose(ss_i.apply_maneuver(man).ecc, 0 * u.one, atol=1e-12 * u.one)
    assert_quantity_allclose(man.get_total_cost(), expected_dv, rtol=1e-5)
    assert_quantity_allclose(man.get_total_time(), expected_total_time, rtol=1e-6)


def test_apply_maneuver_correct_dimensions():
    orb = Orbit.from_vectors(
        Moon,
        [-22681.58976181, 942.47776988, 0] * u.km,
        [-0.04578917, -0.19408599, 0.0] * u.km / u.s,
        Time("2023-08-30 23:14", scale="tdb"),
    )
    man = Maneuver((1 * u.s, [0.01, 0, 0] * u.km / u.s))

    new_orb = orb.apply_maneuver(man, intermediate=False)

    assert new_orb.r.ndim == 1
    assert new_orb.v.ndim == 1


def test_repr_maneuver():
    alt_f = 35781.34857 * u.km
    r = [-6045, -3490, 2500] * u.km
    v = [-3.457, 6.618, 2.533] * u.km / u.s
    alt_b = 503873.0 * u.km
    alt_fi = 376310.0 * u.km
    ss_i = Orbit.from_vectors(Earth, r, v)

    expected_hohmann_maneuver = "Number of impulses: 2, Total cost: 3.060548 km / s"
    expected_bielliptic_maneuver = "Number of impulses: 3, Total cost: 3.122556 km / s"

    assert repr(Maneuver.hohmann(ss_i, Earth.R + alt_f)) == expected_hohmann_maneuver
    assert (
        repr(Maneuver.bielliptic(ss_i, Earth.R + alt_b, Earth.R + alt_fi))
        == expected_bielliptic_maneuver
    )


# Similar Example obtained from "Fundamentals of Astrodynamics and Applications, 4th ed (2013)" by David A. Vallado, page 895
@pytest.mark.parametrize(
    "attractor, max_delta_r, a, ecc, inc, expected_t, expected_v",
    [
        (
            Earth,
            30 * u.km,
            6570 * u.km,
            0.001 * u.one,
            0.7855682278773197 * u.rad,
            2224141.03634 * u.s,
            np.array([0, 0.0083290328315531, 0.00833186625871848]) * (u.km / u.s),
        ),
    ],
)
def test_correct_pericenter(
    attractor, max_delta_r, a, ecc, inc, expected_t, expected_v
):
    ss0 = Orbit.from_classical(
        attractor,
        a,
        ecc,
        inc,
        0 * u.deg,
        0 * u.deg,
        0 * u.deg,
    )

    maneuver = Maneuver.correct_pericenter(ss0, max_delta_r)
    assert_quantity_allclose(maneuver[0][0], expected_t)
    assert_quantity_allclose(maneuver[0][1].value.tolist(), expected_v.value.tolist())


def test_correct_pericenter_J2_exception():
    ss0 = Orbit.from_classical(
        Mercury,
        1000 * u.km,
        0 * u.one,
        0 * u.deg,
        0 * u.deg,
        0 * u.deg,
        0 * u.deg,
    )
    max_delta_r = 30 * u.km
    with pytest.raises(NotImplementedError) as excinfo:
        Maneuver.correct_pericenter(ss0, max_delta_r)
    assert excinfo.type == NotImplementedError
    assert (
        str(excinfo.value)
        == f"The correction maneuver is not yet supported for {ss0.attractor}"
    )


def test_correct_pericenter_ecc_exception():
    ss0 = Orbit.from_classical(
        Earth,
        1000 * u.km,
        0.5 * u.one,
        0 * u.deg,
        0 * u.deg,
        0 * u.deg,
        0 * u.deg,
    )
    max_delta_r = 30 * u.km
    with pytest.raises(NotImplementedError) as excinfo:
        Maneuver.correct_pericenter(ss0, max_delta_r)
    assert excinfo.type == NotImplementedError
    assert (
        str(excinfo.value)
        == f"The correction maneuver is not yet supported with {ss0.ecc},it should be less than or equal to 0.001"
    )
