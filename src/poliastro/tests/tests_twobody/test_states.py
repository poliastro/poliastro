from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from numpy.testing import assert_allclose

from poliastro.bodies import Earth, Sun
from poliastro.twobody._states import ClassicalState, RVState


def test_state_has_attractor_given_in_constructor():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    ss = ClassicalState(Sun, _d, _, _a, _a, _a, _a)
    assert ss.attractor == Sun


def test_state_has_elements_given_in_constructor():
    # Mars data from HORIZONS at J2000
    a = 1.523679 * u.AU
    ecc = 0.093315 * u.one
    inc = 1.85 * u.deg
    raan = 49.562 * u.deg
    argp = 286.537 * u.deg
    nu = 23.33 * u.deg
    ss = ClassicalState(Sun, a * (1 - ecc ** 2), ecc, inc, raan, argp, nu)
    assert ss.coe() == (a, ecc, inc, raan, argp, nu)


def test_state_has_individual_elements():
    a = 1.523679 * u.AU
    ecc = 0.093315 * u.one
    inc = 1.85 * u.deg
    raan = 49.562 * u.deg
    argp = 286.537 * u.deg
    nu = 23.33 * u.deg
    ss = ClassicalState(Sun, a * (1 - ecc ** 2), ecc, inc, raan, argp, nu)
    assert ss.a == a
    assert ss.ecc == ecc
    assert ss.inc == inc
    assert ss.raan == raan
    assert ss.argp == argp
    assert ss.nu == nu


def test_state_has_rv_given_in_constructor():
    r = [1.0, 0.0, 0.0] * u.AU
    v = [0.0, 1.0e-6, 0.0] * u.AU / u.s
    ss = RVState(Sun, r, v)
    assert ss.rv() == (r, v)


def test_perigee_and_apogee():
    expected_r_a = 500 * u.km
    expected_r_p = 300 * u.km
    a = (expected_r_a + expected_r_p) / 2
    ecc = expected_r_a / a - 1
    _a = 1.0 * u.deg  # Unused angle
    ss = ClassicalState(Earth, a * (1 - ecc ** 2), ecc, _a, _a, _a, _a)
    assert_allclose(ss.r_a.to(u.km).value, expected_r_a.to(u.km).value)
    assert_allclose(ss.r_p.to(u.km).value, expected_r_p.to(u.km).value)


def test_convert_from_rv_to_coe():
    # Data from Vallado, example 2.6
    attractor = Earth
    p = 11067.790 * u.km
    ecc = 0.83285 * u.one
    inc = 87.87 * u.deg
    raan = 227.89 * u.deg
    argp = 53.38 * u.deg
    nu = 92.335 * u.deg
    expected_r = [6525.344, 6861.535, 6449.125] * u.km
    expected_v = [4.902276, 5.533124, -1.975709] * u.km / u.s

    r, v = ClassicalState(attractor, p, ecc, inc, raan, argp, nu).rv()

    assert_quantity_allclose(r, expected_r, rtol=1e-5)
    assert_quantity_allclose(v, expected_v, rtol=1e-5)


def test_convert_from_coe_to_rv():
    # Data from Vallado, example 2.5
    attractor = Earth
    r = [6524.384, 6862.875, 6448.296] * u.km
    v = [4.901327, 5.533756, -1.976341] * u.km / u.s

    expected_p = 11067.79 * u.km
    expected_ecc = 0.832853 * u.one
    expected_inc = 87.870 * u.deg
    expected_raan = 227.89 * u.deg
    expected_argp = 53.38 * u.deg
    expected_nu = 92.335 * u.deg

    ss = RVState(attractor, r, v)

    _, ecc, inc, raan, argp, nu = ss.coe()
    p = ss.p

    assert_quantity_allclose(p, expected_p, rtol=1e-4)
    assert_quantity_allclose(ecc, expected_ecc, rtol=1e-4)
    assert_quantity_allclose(inc, expected_inc, rtol=1e-4)
    assert_quantity_allclose(raan, expected_raan, rtol=1e-4)
    assert_quantity_allclose(argp, expected_argp, rtol=1e-4)
    assert_quantity_allclose(nu, expected_nu, rtol=1e-4)


def test_perifocal_points_to_perigee():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    ss = ClassicalState(Sun, _d, _, _a, _a, _a, _a)
    p, _, _ = ss.pqw()
    assert_allclose(p, ss.e_vec / ss.ecc)


def test_arglat_within_range():
    r = [3539.08827417, 5310.19903462, 3066.31301457] * u.km
    v = [-6.49780849, 3.24910291, 1.87521413] * u.km / u.s
    ss = RVState(Earth, r, v)
    assert 0 * u.deg <= ss.arglat <= 360 * u.deg


def test_pqw_returns_dimensionless():
    r_0 = ([1, 0, 0] * u.au).to(u.km)
    v_0 = ([0, 6, 0] * u.au / u.year).to(u.km / u.day)
    ss = RVState(Sun, r_0, v_0)

    p, q, w = ss.pqw()

    assert p.unit == u.one
    assert q.unit == u.one
    assert w.unit == u.one
