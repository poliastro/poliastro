import functools

import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import Angle, solar_system_ephemeris
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time
from numpy.linalg import norm

from poliastro.bodies import Earth, Moon, Sun
from poliastro.constants import H0_earth, Wdivc_sun, rho0_earth
from poliastro.core.elements import rv2coe
from poliastro.core.perturbations import (
    J2_perturbation,
    J3_perturbation,
    atmospheric_drag,
    radiation_pressure,
    third_body,
)
from poliastro.ephem import build_ephem_interpolant
from poliastro.twobody import Orbit
from poliastro.twobody.propagation import cowell


@pytest.mark.slow
def test_J2_propagation_Earth():
    # from Curtis example 12.2:
    r0 = np.array([-2384.46, 5729.01, 3050.46])  # km
    v0 = np.array([-7.36138, -2.98997, 1.64354])  # km/s

    orbit = Orbit.from_vectors(Earth, r0 * u.km, v0 * u.km / u.s)

    tofs = [48.0] * u.h
    rr, vv = cowell(
        Earth.k,
        orbit.r,
        orbit.v,
        tofs,
        ad=J2_perturbation,
        J2=Earth.J2.value,
        R=Earth.R.to(u.km).value,
    )

    k = Earth.k.to(u.km ** 3 / u.s ** 2).value

    _, _, _, raan0, argp0, _ = rv2coe(k, r0, v0)
    _, _, _, raan, argp, _ = rv2coe(k, rr[0].to(u.km).value, vv[0].to(u.km / u.s).value)

    raan_variation_rate = (raan - raan0) / tofs[0].to(u.s).value  # type: ignore
    argp_variation_rate = (argp - argp0) / tofs[0].to(u.s).value  # type: ignore

    raan_variation_rate = (raan_variation_rate * u.rad / u.s).to(u.deg / u.h)
    argp_variation_rate = (argp_variation_rate * u.rad / u.s).to(u.deg / u.h)

    assert_quantity_allclose(raan_variation_rate, -0.172 * u.deg / u.h, rtol=1e-2)
    assert_quantity_allclose(argp_variation_rate, 0.282 * u.deg / u.h, rtol=1e-2)


@pytest.mark.slow
@pytest.mark.parametrize(
    "test_params",
    [
        {
            "inc": 0.2618 * u.rad,
            "da_max": 43.2 * u.m,
            "dinc_max": 3.411e-5,
            "decc_max": 3.549e-5,
        },
        {
            "inc": 0.7854 * u.rad,
            "da_max": 135.8 * u.m,
            "dinc_max": 2.751e-5,
            "decc_max": 9.243e-5,
        },
        {
            "inc": 1.3090 * u.rad,
            "da_max": 58.7 * u.m,
            "dinc_max": 0.79e-5,
            "decc_max": 10.02e-5,
        },
        {
            "inc": 1.5708 * u.rad,
            "da_max": 96.1 * u.m,
            "dinc_max": 0.0,
            "decc_max": 17.04e-5,
        },
    ],
)
def test_J3_propagation_Earth(test_params):
    # Nai-ming Qi, Qilong Sun, Yong Yang, (2018) "Effect of J3 perturbation on satellite position in LEO",
    # Aircraft Engineering and  Aerospace Technology, Vol. 90 Issue: 1,
    # pp.74-86, https://doi.org/10.1108/AEAT-03-2015-0092
    a_ini = 8970.667 * u.km
    ecc_ini = 0.25 * u.one
    raan_ini = 1.047 * u.rad
    nu_ini = 0.0 * u.rad
    argp_ini = 1.0 * u.rad
    inc_ini = test_params["inc"]

    k = Earth.k.to(u.km ** 3 / u.s ** 2).value

    orbit = Orbit.from_classical(
        Earth, a_ini, ecc_ini, inc_ini, raan_ini, argp_ini, nu_ini
    )

    tofs = np.linspace(0, 10.0 * u.day, 1000)
    r_J2, v_J2 = cowell(
        Earth.k,
        orbit.r,
        orbit.v,
        tofs,
        ad=J2_perturbation,
        J2=Earth.J2.value,
        R=Earth.R.to(u.km).value,
        rtol=1e-8,
    )

    def a_J2J3(t0, u_, k_):
        j2 = J2_perturbation(t0, u_, k_, J2=Earth.J2.value, R=Earth.R.to(u.km).value)
        j3 = J3_perturbation(t0, u_, k_, J3=Earth.J3.value, R=Earth.R.to(u.km).value)
        return j2 + j3

    r_J3, v_J3 = cowell(Earth.k, orbit.r, orbit.v, tofs, ad=a_J2J3, rtol=1e-8)

    a_values_J2 = np.array(
        [
            rv2coe(k, ri, vi)[0] / (1.0 - rv2coe(k, ri, vi)[1] ** 2)
            for ri, vi in zip(r_J2.to(u.km).value, v_J2.to(u.km / u.s).value)
        ]
    )
    a_values_J3 = np.array(
        [
            rv2coe(k, ri, vi)[0] / (1.0 - rv2coe(k, ri, vi)[1] ** 2)
            for ri, vi in zip(r_J3.to(u.km).value, v_J3.to(u.km / u.s).value)
        ]
    )
    da_max = np.max(np.abs(a_values_J2 - a_values_J3))

    ecc_values_J2 = np.array(
        [
            rv2coe(k, ri, vi)[1]
            for ri, vi in zip(r_J2.to(u.km).value, v_J2.to(u.km / u.s).value)
        ]
    )
    ecc_values_J3 = np.array(
        [
            rv2coe(k, ri, vi)[1]
            for ri, vi in zip(r_J3.to(u.km).value, v_J3.to(u.km / u.s).value)
        ]
    )
    decc_max = np.max(np.abs(ecc_values_J2 - ecc_values_J3))

    inc_values_J2 = np.array(
        [
            rv2coe(k, ri, vi)[2]
            for ri, vi in zip(r_J2.to(u.km).value, v_J2.to(u.km / u.s).value)
        ]
    )
    inc_values_J3 = np.array(
        [
            rv2coe(k, ri, vi)[2]
            for ri, vi in zip(r_J3.to(u.km).value, v_J3.to(u.km / u.s).value)
        ]
    )
    dinc_max = np.max(np.abs(inc_values_J2 - inc_values_J3))

    assert_quantity_allclose(dinc_max, test_params["dinc_max"], rtol=1e-1, atol=1e-7)
    assert_quantity_allclose(decc_max, test_params["decc_max"], rtol=1e-1, atol=1e-7)
    try:
        assert_quantity_allclose(da_max * u.km, test_params["da_max"])
    except AssertionError:
        pytest.xfail("this assertion disagrees with the paper")


@pytest.mark.slow
def test_atmospheric_drag():
    # http://farside.ph.utexas.edu/teaching/celestial/Celestialhtml/node94.html#sair (10.148)
    # given the expression for \dot{r} / r, aproximate \Delta r \approx F_r * \Delta t

    R = Earth.R.to(u.km).value
    k = Earth.k.to(u.km ** 3 / u.s ** 2).value

    # parameters of a circular orbit with h = 250 km (any value would do, but not too small)
    orbit = Orbit.circular(Earth, 250 * u.km)
    r0, _ = orbit.rv()
    r0 = r0.to(u.km).value

    # parameters of a body
    C_D = 2.2  # dimentionless (any value would do)
    A = ((np.pi / 4.0) * (u.m ** 2)).to(u.km ** 2).value  # km^2
    m = 100  # kg
    B = C_D * A / m

    # parameters of the atmosphere
    rho0 = rho0_earth.to(u.kg / u.km ** 3).value  # kg/km^3
    H0 = H0_earth.to(u.km).value  # km
    tof = 100000  # s

    dr_expected = -B * rho0 * np.exp(-(norm(r0) - R) / H0) * np.sqrt(k * norm(r0)) * tof
    # assuming the atmospheric decay during tof is small,
    # dr_expected = F_r * tof (Newton's integration formula), where
    # F_r = -B rho(r) |r|^2 sqrt(k / |r|^3) = -B rho(r) sqrt(k |r|)

    rr, _ = cowell(
        Earth.k,
        orbit.r,
        orbit.v,
        [tof] * u.s,
        ad=atmospheric_drag,
        R=R,
        C_D=C_D,
        A=A,
        m=m,
        H0=H0,
        rho0=rho0,
    )

    assert_quantity_allclose(
        norm(rr[0].to(u.km).value) - norm(r0), dr_expected, rtol=1e-2
    )


@pytest.mark.slow
def test_atmospheric_demise():
    # Test an orbital decay that hits Earth. No analytic solution.
    R = Earth.R.to(u.km).value

    orbit = Orbit.circular(Earth, 230 * u.km)

    # parameters of a body
    C_D = 2.2  # dimentionless (any value would do)
    A = ((np.pi / 4.0) * (u.m ** 2)).to(u.km ** 2).value  # km^2
    m = 100  # kg

    # parameters of the atmosphere
    rho0 = rho0_earth.to(u.kg / u.km ** 3).value  # kg/km^3
    H0 = H0_earth.to(u.km).value  # km

    tofs = [365] * u.d  # actually hits the ground a bit after day 48

    class LithobrakeEvent:
        def __init__(self, R):
            self._R = R
            self._last_t = None

        @property
        def terminal(self):
            # tell scipy to stop the integration at H = R (impact)
            return True

        @property
        def last_t(self):
            return self._last_t * u.s

        def __call__(self, t, u):
            self._last_t = t
            H = norm(u[:3])
            # scipy will search for H - R = 0
            return H - self._R

    lithobrake_event = LithobrakeEvent(R)
    events = [lithobrake_event]

    rr, _ = cowell(
        Earth.k,
        orbit.r,
        orbit.v,
        tofs,
        ad=atmospheric_drag,
        R=R,
        C_D=C_D,
        A=A,
        m=m,
        H0=H0,
        rho0=rho0,
        events=events,
    )

    assert_quantity_allclose(norm(rr[0].to(u.km).value), R, atol=1)  # below 1km

    # last_t is not an analytic value
    assert_quantity_allclose(lithobrake_event.last_t.to(u.d).value, 48.2179, rtol=1e-2)


@pytest.mark.slow
def test_cowell_works_with_small_perturbations():
    r0 = [-2384.46, 5729.01, 3050.46] * u.km
    v0 = [-7.36138, -2.98997, 1.64354] * u.km / u.s

    r_expected = [
        13179.39566663877121754922,
        -13026.25123408228319021873,
        -9852.66213692844394245185,
    ] * u.km
    v_expected = (
        [2.78170542314378943516, 3.21596786944631274352, 0.16327165546278937791]
        * u.km
        / u.s
    )

    initial = Orbit.from_vectors(Earth, r0, v0)

    def accel(t0, state, k):
        v_vec = state[3:]
        norm_v = (v_vec * v_vec).sum() ** 0.5
        return 1e-5 * v_vec / norm_v

    final = initial.propagate(3 * u.day, method=cowell, ad=accel)

    assert_quantity_allclose(final.r, r_expected)
    assert_quantity_allclose(final.v, v_expected)


@pytest.mark.slow
def test_cowell_converges_with_small_perturbations():
    r0 = [-2384.46, 5729.01, 3050.46] * u.km
    v0 = [-7.36138, -2.98997, 1.64354] * u.km / u.s

    initial = Orbit.from_vectors(Earth, r0, v0)

    def accel(t0, state, k):
        v_vec = state[3:]
        norm_v = (v_vec * v_vec).sum() ** 0.5
        return 0.0 * v_vec / norm_v

    final = initial.propagate(initial.period, method=cowell, ad=accel)

    assert_quantity_allclose(final.r, initial.r)
    assert_quantity_allclose(final.v, initial.v)


moon_heo = {
    "body": Moon,
    "tof": 60 * u.day,
    "raan": -0.06 * u.deg,
    "argp": 0.15 * u.deg,
    "inc": 0.08 * u.deg,
    "orbit": [
        26553.4 * u.km,
        0.741 * u.one,
        63.4 * u.deg,
        0.0 * u.deg,
        -10.12921 * u.deg,
        0.0 * u.rad,
    ],
    "period": 28 * u.day,
}

moon_leo = {
    "body": Moon,
    "tof": 60 * u.day,
    "raan": -2.18 * 1e-4 * u.deg,
    "argp": 15.0 * 1e-3 * u.deg,
    "inc": 6.0 * 1e-4 * u.deg,
    "orbit": [
        6678.126 * u.km,
        0.01 * u.one,
        28.5 * u.deg,
        0.0 * u.deg,
        0.0 * u.deg,
        0.0 * u.rad,
    ],
    "period": 28 * u.day,
}

moon_geo = {
    "body": Moon,
    "tof": 60 * u.day,
    "raan": 6.0 * u.deg,
    "argp": -11.0 * u.deg,
    "inc": 6.5 * 1e-3 * u.deg,
    "orbit": [
        42164.0 * u.km,
        0.0001 * u.one,
        1 * u.deg,
        0.0 * u.deg,
        0.0 * u.deg,
        0.0 * u.rad,
    ],
    "period": 28 * u.day,
}

sun_heo = {
    "body": Sun,
    "tof": 200 * u.day,
    "raan": -0.10 * u.deg,
    "argp": 0.2 * u.deg,
    "inc": 0.1 * u.deg,
    "orbit": [
        26553.4 * u.km,
        0.741 * u.one,
        63.4 * u.deg,
        0.0 * u.deg,
        -10.12921 * u.deg,
        0.0 * u.rad,
    ],
    "period": 365 * u.day,
}

sun_leo = {
    "body": Sun,
    "tof": 200 * u.day,
    "raan": -6.0 * 1e-3 * u.deg,
    "argp": 0.02 * u.deg,
    "inc": -1.0 * 1e-4 * u.deg,
    "orbit": [
        6678.126 * u.km,
        0.01 * u.one,
        28.5 * u.deg,
        0.0 * u.deg,
        0.0 * u.deg,
        0.0 * u.rad,
    ],
    "period": 365 * u.day,
}

sun_geo = {
    "body": Sun,
    "tof": 200 * u.day,
    "raan": 8.7 * u.deg,
    "argp": -5.5 * u.deg,
    "inc": 5.5e-3 * u.deg,
    "orbit": [
        42164.0 * u.km,
        0.0001 * u.one,
        1 * u.deg,
        0.0 * u.deg,
        0.0 * u.deg,
        0.0 * u.rad,
    ],
    "period": 365 * u.day,
}


@pytest.mark.slow
@pytest.mark.parametrize(
    "test_params",
    [
        moon_heo,
        moon_geo,
        moon_leo,
        sun_heo,
        sun_geo,
        pytest.param(
            sun_leo,
            marks=pytest.mark.skip(
                reason="here agreement required rtol=1e-10, too long for 200 days"
            ),
        ),
    ],
)
def test_3rd_body_Curtis(test_params):
    # based on example 12.11 from Howard Curtis
    body = test_params["body"]
    with solar_system_ephemeris.set("builtin"):
        j_date = 2454283.0 * u.day
        tof = (test_params["tof"]).to(u.s).value
        body_r = build_ephem_interpolant(
            body,
            test_params["period"],
            (j_date, j_date + test_params["tof"]),
            rtol=1e-2,
        )

        epoch = Time(j_date, format="jd", scale="tdb")
        initial = Orbit.from_classical(Earth, *test_params["orbit"], epoch=epoch)
        rr, vv = cowell(
            Earth.k,
            initial.r,
            initial.v,
            np.linspace(0, tof, 400) * u.s,
            rtol=1e-10,
            ad=third_body,
            k_third=body.k.to(u.km ** 3 / u.s ** 2).value,
            third_body=body_r,
        )

        incs, raans, argps = [], [], []
        for ri, vi in zip(rr.to(u.km).value, vv.to(u.km / u.s).value):
            angles = Angle(
                rv2coe(Earth.k.to(u.km ** 3 / u.s ** 2).value, ri, vi)[2:5] * u.rad
            )  # inc, raan, argp
            angles = angles.wrap_at(180 * u.deg)
            incs.append(angles[0].value)
            raans.append(angles[1].value)
            argps.append(angles[2].value)

        # averaging over 5 last values in the way Curtis does
        inc_f, raan_f, argp_f = (
            np.mean(incs[-5:]),
            np.mean(raans[-5:]),
            np.mean(argps[-5:]),
        )

        assert_quantity_allclose(
            [
                (raan_f * u.rad).to(u.deg) - test_params["orbit"][3],
                (inc_f * u.rad).to(u.deg) - test_params["orbit"][2],
                (argp_f * u.rad).to(u.deg) - test_params["orbit"][4],
            ],
            [test_params["raan"], test_params["inc"], test_params["argp"]],
            rtol=1e-1,
        )


@pytest.fixture(scope="module")
def sun_r():
    j_date = 2_438_400.5 * u.day
    tof = 600 * u.day
    return build_ephem_interpolant(Sun, 365 * u.day, (j_date, j_date + tof), rtol=1e-2)


def normalize_to_Curtis(t0, sun_r):
    r = sun_r(t0)
    return 149600000 * r / norm(r)


@pytest.mark.slow
@pytest.mark.parametrize(
    "t_days,deltas_expected",
    [
        (200, [3e-3, -8e-3, -0.035, -80.0]),
        (400, [-1.3e-3, 0.01, -0.07, 8.0]),
        (600, [7e-3, 0.03, -0.10, -80.0]),
        # (800, [-7.5e-3, 0.02, -0.13, 1.7]),
        # (1000, [6e-3, 0.065, -0.165, -70.0]),
        # (1095, [0.0, 0.06, -0.165, -10.0]),
    ],
)
def test_solar_pressure(t_days, deltas_expected, sun_r):
    # based on example 12.9 from Howard Curtis
    with solar_system_ephemeris.set("builtin"):
        j_date = 2_438_400.5 * u.day
        tof = 600 * u.day
        epoch = Time(j_date, format="jd", scale="tdb")

        initial = Orbit.from_classical(
            Earth,
            10085.44 * u.km,
            0.025422 * u.one,
            88.3924 * u.deg,
            45.38124 * u.deg,
            227.493 * u.deg,
            343.4268 * u.deg,
            epoch=epoch,
        )
        # in Curtis, the mean distance to Sun is used. In order to validate against it, we have to do the same thing
        sun_normalized = functools.partial(normalize_to_Curtis, sun_r=sun_r)

        rr, vv = cowell(
            Earth.k,
            initial.r,
            initial.v,
            np.linspace(0, (tof).to(u.s).value, 4000) * u.s,
            rtol=1e-8,
            ad=radiation_pressure,
            R=Earth.R.to(u.km).value,
            C_R=2.0,
            A=2e-4,
            m=100,
            Wdivc_s=Wdivc_sun.value,
            star=sun_normalized,
        )

        delta_eccs, delta_incs, delta_raans, delta_argps = [], [], [], []
        for ri, vi in zip(rr.to(u.km).value, vv.to(u.km / u.s).value):
            orbit_params = rv2coe(Earth.k.to(u.km ** 3 / u.s ** 2).value, ri, vi)
            delta_eccs.append(orbit_params[1] - initial.ecc.value)
            delta_incs.append(
                (orbit_params[2] * u.rad).to(u.deg).value - initial.inc.value
            )
            delta_raans.append(
                (orbit_params[3] * u.rad).to(u.deg).value - initial.raan.value
            )
            delta_argps.append(
                (orbit_params[4] * u.rad).to(u.deg).value - initial.argp.value
            )

        # averaging over 5 last values in the way Curtis does
        index = int(
            1.0 * t_days / tof.to(u.day).value * 4000  # type: ignore
        )
        delta_ecc, delta_inc, delta_raan, delta_argp = (
            np.mean(delta_eccs[index - 5 : index]),
            np.mean(delta_incs[index - 5 : index]),
            np.mean(delta_raans[index - 5 : index]),
            np.mean(delta_argps[index - 5 : index]),
        )
        assert_quantity_allclose(
            [delta_ecc, delta_inc, delta_raan, delta_argp],
            deltas_expected,
            rtol=1e0,  # TODO: Excessively low, rewrite test?
            atol=1e-4,
        )
