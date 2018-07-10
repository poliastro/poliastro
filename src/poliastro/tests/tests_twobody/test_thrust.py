import pytest

import numpy as np
from numpy.testing import assert_allclose

from astropy import units as u

from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.twobody.propagation import cowell
from poliastro.twobody.thrust import change_a_inc, change_argp


def test_soyuz_standard_gto_delta_v():
    # Data from Soyuz Users Manual, issue 2 revision 0
    r_a = (Earth.R + 35950 * u.km).to(u.km).value
    r_p = (Earth.R + 250 * u.km).to(u.km).value

    a = (r_a + r_p) / 2  # km
    ecc = r_a / a - 1
    argp_0 = (178 * u.deg).to(u.rad).value  # rad
    argp_f = (178 * u.deg + 5 * u.deg).to(u.rad).value  # rad
    f = 2.4e-7  # km / s2

    k = Earth.k.to(u.km**3 / u.s**2).value

    _, delta_V, t_f = change_argp(k, a, ecc, argp_0, argp_f, f)

    expected_t_f = 12.0  # days, approximate
    expected_delta_V = 0.2489  # km / s

    assert_allclose(delta_V, expected_delta_V, rtol=1e-2)
    assert_allclose(t_f / 86400, expected_t_f, rtol=1e-2)


def test_soyuz_standard_gto_numerical():
    # Data from Soyuz Users Manual, issue 2 revision 0
    r_a = (Earth.R + 35950 * u.km).to(u.km).value
    r_p = (Earth.R + 250 * u.km).to(u.km).value

    a = (r_a + r_p) / 2  # km
    ecc = r_a / a - 1
    argp_0 = (178 * u.deg).to(u.rad).value  # rad
    argp_f = (178 * u.deg + 5 * u.deg).to(u.rad).value  # rad
    f = 2.4e-7  # km / s2

    k = Earth.k.to(u.km**3 / u.s**2).value

    argp_accel, _, t_f = change_argp(k, a, ecc, argp_0, argp_f, f)

    # Retrieve r and v from initial orbit
    s0 = Orbit.from_classical(
        Earth,
        a * u.km, (r_a / a - 1) * u.one, 6 * u.deg,
        188.5 * u.deg, 178 * u.deg, 0 * u.deg
    )

    # Propagate orbit
    r, v = cowell(s0, t_f, ad=argp_accel, rtol=1e-8)

    sf = Orbit.from_vectors(Earth,
                            r * u.km,
                            v * u.km / u.s,
                            s0.epoch + t_f * u.s)

    assert_allclose(sf.argp.to(u.rad).value, argp_f, rtol=1e-4)


@pytest.mark.parametrize("inc_0, expected_t_f, expected_delta_V, rtol", [
    [28.5, 191.26295, 5.78378, 1e-5],
    [90.0, 335.0, 10.13, 1e-3],
    [114.591, 351.0, 10.61, 1e-2]
])
def test_leo_geo_time_and_delta_v(inc_0, expected_t_f, expected_delta_V, rtol):
    f = 3.5e-7  # km / s2

    a_0 = 7000.0  # km
    a_f = 42166.0  # km
    inc_f = 0.0  # rad
    k = Earth.k.to(u.km**3 / u.s**2).value
    inc_0 = np.radians(inc_0)  # rad

    _, delta_V, t_f = change_a_inc(k, a_0, a_f, inc_0, inc_f, f)

    assert_allclose(delta_V, expected_delta_V, rtol=rtol)
    assert_allclose((t_f * u.s).to(u.day).value, expected_t_f, rtol=rtol)


@pytest.mark.parametrize("inc_0", [np.radians(28.5), np.radians(90.0)])
def test_leo_geo_numerical(inc_0):
    f = 3.5e-7  # km / s2

    a_0 = 7000.0  # km
    a_f = 42166.0  # km
    inc_f = 0.0  # rad

    k = Earth.k.to(u.km**3 / u.s**2).value

    edelbaum_accel, _, t_f = change_a_inc(k, a_0, a_f, inc_0, inc_f, f)

    # Retrieve r and v from initial orbit
    s0 = Orbit.circular(Earth, a_0 * u.km - Earth.R, inc_0 * u.rad)
    r0, v0 = s0.rv()

    # Propagate orbit

    r, v = cowell(s0, t_f, ad=edelbaum_accel, rtol=1e-7)

    sf = Orbit.from_vectors(Earth, r * u.km, v * u.km / u.s, s0.epoch + t_f * u.s)

    assert_allclose(sf.a.to(u.km).value, a_f, rtol=1e-3)
    assert_allclose(sf.ecc.value, 0.0, atol=1e-2)
    assert_allclose(sf.inc.to(u.rad).value, inc_f, atol=2e-3)
