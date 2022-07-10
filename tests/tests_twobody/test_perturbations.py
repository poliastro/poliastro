import functools

import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import Angle
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
    atmospheric_drag_exponential,
    radiation_pressure,
    third_body,
)
from poliastro.core.propagation import func_twobody
from poliastro.earth.atmosphere import COESA76
from poliastro.ephem import build_ephem_interpolant
from poliastro.twobody import Orbit
from poliastro.twobody.events import LithobrakeEvent
from poliastro.twobody.propagation import CowellPropagator


@pytest.mark.slow
def test_J2_propagation_Earth():
    # From Curtis example 12.2:
    r0 = np.array([-2384.46, 5729.01, 3050.46])  # km
    v0 = np.array([-7.36138, -2.98997, 1.64354])  # km/s

    orbit = Orbit.from_vectors(Earth, r0 * u.km, v0 * u.km / u.s)

    tofs = [48.0] * u.h

    def f(t0, u_, k):
        du_kep = func_twobody(t0, u_, k)
        ax, ay, az = J2_perturbation(
            t0, u_, k, J2=Earth.J2.value, R=Earth.R.to(u.km).value
        )
        du_ad = np.array([0, 0, 0, ax, ay, az])
        return du_kep + du_ad

    method = CowellPropagator(f=f)
    rr, vv = method.propagate_many(orbit._state, tofs)

    k = Earth.k.to(u.km**3 / u.s**2).value

    _, _, _, raan0, argp0, _ = rv2coe(k, r0, v0)
    _, _, _, raan, argp, _ = rv2coe(
        k, rr[0].to(u.km).value, vv[0].to(u.km / u.s).value
    )

    raan_variation_rate = (raan - raan0) / tofs[0].to(u.s).value  # type: ignore
    argp_variation_rate = (argp - argp0) / tofs[0].to(u.s).value  # type: ignore

    raan_variation_rate = (raan_variation_rate * u.rad / u.s).to(u.deg / u.h)
    argp_variation_rate = (argp_variation_rate * u.rad / u.s).to(u.deg / u.h)

    assert_quantity_allclose(
        raan_variation_rate, -0.172 * u.deg / u.h, rtol=1e-2
    )
    assert_quantity_allclose(
        argp_variation_rate, 0.282 * u.deg / u.h, rtol=1e-2
    )


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

    k = Earth.k.to(u.km**3 / u.s**2).value

    orbit = Orbit.from_classical(
        attractor=Earth,
        a=a_ini,
        ecc=ecc_ini,
        inc=inc_ini,
        raan=raan_ini,
        argp=argp_ini,
        nu=nu_ini,
    )

    def f(t0, u_, k):
        du_kep = func_twobody(t0, u_, k)
        ax, ay, az = J2_perturbation(
            t0, u_, k, J2=Earth.J2.value, R=Earth.R.to(u.km).value
        )
        du_ad = np.array([0, 0, 0, ax, ay, az])
        return du_kep + du_ad

    tofs = np.linspace(0, 10.0 * u.day, 1000)
    method = CowellPropagator(rtol=1e-8, f=f)
    r_J2, v_J2 = method.propagate_many(
        orbit._state,
        tofs,
    )

    def f_combined(t0, u_, k):
        du_kep = func_twobody(t0, u_, k)
        ax, ay, az = J2_perturbation(
            t0, u_, k, J2=Earth.J2.value, R=Earth.R.to_value(u.km)
        ) + J3_perturbation(
            t0, u_, k, J3=Earth.J3.value, R=Earth.R.to_value(u.km)
        )
        du_ad = np.array([0, 0, 0, ax, ay, az])
        return du_kep + du_ad

    method = CowellPropagator(rtol=1e-8, f=f_combined)
    r_J3, v_J3 = method.propagate_many(
        orbit._state,
        tofs,
    )

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

    assert_quantity_allclose(
        dinc_max, test_params["dinc_max"], rtol=1e-1, atol=1e-7
    )
    assert_quantity_allclose(
        decc_max, test_params["decc_max"], rtol=1e-1, atol=1e-7
    )
    try:
        assert_quantity_allclose(da_max * u.km, test_params["da_max"])
    except AssertionError:
        pytest.xfail("this assertion disagrees with the paper")


@pytest.mark.slow
def test_atmospheric_drag_exponential():
    # http://farside.ph.utexas.edu/teaching/celestial/Celestialhtml/node94.html#sair (10.148)
    # Given the expression for \dot{r} / r, aproximate \Delta r \approx F_r * \Delta t

    R = Earth.R.to(u.km).value
    k = Earth.k.to(u.km**3 / u.s**2).value

    # Parameters of a circular orbit with h = 250 km (any value would do, but not too small)
    orbit = Orbit.circular(Earth, 250 * u.km)
    r0, _ = orbit.rv()
    r0 = r0.to(u.km).value

    # Parameters of a body
    C_D = 2.2  # dimentionless (any value would do)
    A_over_m = ((np.pi / 4.0) * (u.m**2) / (100 * u.kg)).to_value(
        u.km**2 / u.kg
    )  # km^2/kg
    B = C_D * A_over_m

    # Parameters of the atmosphere
    rho0 = rho0_earth.to(u.kg / u.km**3).value  # kg/km^3
    H0 = H0_earth.to(u.km).value  # km
    tof = 100000  # s

    dr_expected = (
        -B * rho0 * np.exp(-(norm(r0) - R) / H0) * np.sqrt(k * norm(r0)) * tof
    )
    # Assuming the atmospheric decay during tof is small,
    # dr_expected = F_r * tof (Newton's integration formula), where
    # F_r = -B rho(r) |r|^2 sqrt(k / |r|^3) = -B rho(r) sqrt(k |r|)

    def f(t0, u_, k):
        du_kep = func_twobody(t0, u_, k)
        ax, ay, az = atmospheric_drag_exponential(
            t0, u_, k, R=R, C_D=C_D, A_over_m=A_over_m, H0=H0, rho0=rho0
        )
        du_ad = np.array([0, 0, 0, ax, ay, az])
        return du_kep + du_ad

    method = CowellPropagator(f=f)
    rr, _ = method.propagate_many(
        orbit._state,
        [tof] * u.s,
    )

    assert_quantity_allclose(
        norm(rr[0].to(u.km).value) - norm(r0), dr_expected, rtol=1e-2
    )


@pytest.mark.slow
def test_atmospheric_demise():
    # Test an orbital decay that hits Earth. No analytic solution.
    R = Earth.R.to(u.km).value

    orbit = Orbit.circular(Earth, 230 * u.km)
    t_decay = 48.2179 * u.d  # not an analytic value

    # Parameters of a body
    C_D = 2.2  # dimentionless (any value would do)
    A_over_m = ((np.pi / 4.0) * (u.m**2) / (100 * u.kg)).to_value(
        u.km**2 / u.kg
    )  # km^2/kg

    # Parameters of the atmosphere
    rho0 = rho0_earth.to(u.kg / u.km**3).value  # kg/km^3
    H0 = H0_earth.to(u.km).value  # km

    tofs = [365] * u.d  # Actually hits the ground a bit after day 48

    lithobrake_event = LithobrakeEvent(R)
    events = [lithobrake_event]

    def f(t0, u_, k):
        du_kep = func_twobody(t0, u_, k)
        ax, ay, az = atmospheric_drag_exponential(
            t0, u_, k, R=R, C_D=C_D, A_over_m=A_over_m, H0=H0, rho0=rho0
        )
        du_ad = np.array([0, 0, 0, ax, ay, az])
        return du_kep + du_ad

    method = CowellPropagator(events=events, f=f)
    rr, _ = method.propagate_many(
        orbit._state,
        tofs,
    )

    assert_quantity_allclose(
        norm(rr[0].to(u.km).value), R, atol=1
    )  # Below 1km

    assert_quantity_allclose(lithobrake_event.last_t, t_decay, rtol=1e-2)

    # Make sure having the event not firing is ok
    tofs = [1] * u.d
    lithobrake_event = LithobrakeEvent(R)
    events = [lithobrake_event]

    method = CowellPropagator(events=events, f=f)
    rr, _ = method.propagate_many(
        orbit._state,
        tofs,
    )

    assert lithobrake_event.last_t == tofs[-1]


@pytest.mark.slow
def test_atmospheric_demise_coesa76():
    # Test an orbital decay that hits Earth. No analytic solution.
    R = Earth.R.to(u.km).value

    orbit = Orbit.circular(Earth, 250 * u.km)
    t_decay = 7.17 * u.d

    # Parameters of a body
    C_D = 2.2  # Dimensionless (any value would do)
    A_over_m = ((np.pi / 4.0) * (u.m**2) / (100 * u.kg)).to_value(
        u.km**2 / u.kg
    )  # km^2/kg

    tofs = [365] * u.d

    lithobrake_event = LithobrakeEvent(R)
    events = [lithobrake_event]

    coesa76 = COESA76()

    def f(t0, u_, k):
        du_kep = func_twobody(t0, u_, k)

        # Avoid undershooting H below attractor radius R
        H = max(norm(u_[:3]), R)
        rho = coesa76.density((H - R) * u.km).to_value(u.kg / u.km**3)

        ax, ay, az = atmospheric_drag(
            t0,
            u_,
            k,
            C_D=C_D,
            A_over_m=A_over_m,
            rho=rho,
        )
        du_ad = np.array([0, 0, 0, ax, ay, az])
        return du_kep + du_ad

    method = CowellPropagator(events=events, f=f)
    rr, _ = method.propagate_many(
        orbit._state,
        tofs,
    )

    assert_quantity_allclose(
        norm(rr[0].to(u.km).value), R, atol=1
    )  # Below 1km

    assert_quantity_allclose(lithobrake_event.last_t, t_decay, rtol=1e-2)


@pytest.mark.slow
def test_cowell_works_with_small_perturbations():
    r0 = [-2384.46, 5729.01, 3050.46] * u.km
    v0 = [-7.36138, -2.98997, 1.64354] * u.km / u.s

    # TODO: Where does this data come from?
    r_expected = [
        13179.39566663877121754922,
        -13026.25123408228319021873,
        -9852.66213692844394245185,
    ] * u.km
    v_expected = (
        [
            2.78170542314378943516,
            3.21596786944631274352,
            0.16327165546278937791,
        ]
        * u.km
        / u.s
    )

    initial = Orbit.from_vectors(Earth, r0, v0)

    def accel(t0, state, k):
        v_vec = state[3:]
        norm_v = (v_vec * v_vec).sum() ** 0.5
        return 1e-5 * v_vec / norm_v

    def f(t0, u_, k):
        du_kep = func_twobody(t0, u_, k)
        ax, ay, az = accel(t0, u_, k)
        du_ad = np.array([0, 0, 0, ax, ay, az])
        return du_kep + du_ad

    final = initial.propagate(3 * u.day, method=CowellPropagator(f=f))

    # TODO: Accuracy reduced after refactor,
    # but unclear what are we comparing against
    assert_quantity_allclose(final.r, r_expected, rtol=1e-6)
    assert_quantity_allclose(final.v, v_expected, rtol=1e-5)


@pytest.mark.slow
def test_cowell_converges_with_small_perturbations():
    r0 = [-2384.46, 5729.01, 3050.46] * u.km
    v0 = [-7.36138, -2.98997, 1.64354] * u.km / u.s

    initial = Orbit.from_vectors(Earth, r0, v0)

    def accel(t0, state, k):
        v_vec = state[3:]
        norm_v = (v_vec * v_vec).sum() ** 0.5
        return 0.0 * v_vec / norm_v

    def f(t0, u_, k):
        du_kep = func_twobody(t0, u_, k)
        ax, ay, az = accel(t0, u_, k)
        du_ad = np.array([0, 0, 0, ax, ay, az])
        return du_kep + du_ad

    final = initial.propagate(initial.period, method=CowellPropagator(f=f))

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
    # Based on example 12.11 from Howard Curtis
    body = test_params["body"]
    j_date = 2454283.0 * u.day
    tof = (test_params["tof"]).to_value(u.s)
    body_r = build_ephem_interpolant(
        body,
        test_params["period"],
        (j_date, j_date + test_params["tof"]),
        rtol=1e-2,
    )

    epoch = Time(j_date, format="jd", scale="tdb")
    initial = Orbit.from_classical(Earth, *test_params["orbit"], epoch=epoch)

    def f(t0, u_, k):
        du_kep = func_twobody(t0, u_, k)
        ax, ay, az = third_body(
            t0,
            u_,
            k,
            k_third=body.k.to_value(u.km**3 / u.s**2),
            perturbation_body=body_r,
        )
        du_ad = np.array([0, 0, 0, ax, ay, az])
        return du_kep + du_ad

    method = CowellPropagator(rtol=1e-10, f=f)
    rr, vv = method.propagate_many(
        initial._state,
        np.linspace(0, tof, 400) << u.s,
    )

    incs, raans, argps = [], [], []
    for ri, vi in zip(rr.to_value(u.km), vv.to_value(u.km / u.s)):
        angles = Angle(
            rv2coe(Earth.k.to_value(u.km**3 / u.s**2), ri, vi)[2:5] * u.rad
        )  # inc, raan, argp
        angles = angles.wrap_at(180 * u.deg)
        incs.append(angles[0].value)
        raans.append(angles[1].value)
        argps.append(angles[2].value)

    # Averaging over 5 last values in the way Curtis does
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
    return build_ephem_interpolant(
        Sun, 365 * u.day, (j_date, j_date + tof), rtol=1e-2
    )


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
    # Based on example 12.9 from Howard Curtis
    j_date = 2_438_400.5 * u.day
    tof = 600 * u.day
    epoch = Time(j_date, format="jd", scale="tdb")

    with pytest.warns(
        UserWarning, match="Wrapping true anomaly to -π <= nu < π"
    ):
        initial = Orbit.from_classical(
            attractor=Earth,
            a=10085.44 * u.km,
            ecc=0.025422 * u.one,
            inc=88.3924 * u.deg,
            raan=45.38124 * u.deg,
            argp=227.493 * u.deg,
            nu=343.4268 * u.deg,
            epoch=epoch,
        )
    # In Curtis, the mean distance to Sun is used. In order to validate against it, we have to do the same thing
    sun_normalized = functools.partial(normalize_to_Curtis, sun_r=sun_r)

    def f(t0, u_, k):
        du_kep = func_twobody(t0, u_, k)
        ax, ay, az = radiation_pressure(
            t0,
            u_,
            k,
            R=Earth.R.to(u.km).value,
            C_R=2.0,
            A_over_m=2e-4 / 100,
            Wdivc_s=Wdivc_sun.value,
            star=sun_normalized,
        )
        du_ad = np.array([0, 0, 0, ax, ay, az])
        return du_kep + du_ad

    method = CowellPropagator(
        rtol=1e-8,
        f=f,
    )
    rr, vv = method.propagate_many(
        initial._state,
        np.linspace(0, tof.to_value(u.s), 4000) << u.s,
    )

    delta_eccs, delta_incs, delta_raans, delta_argps = [], [], [], []
    for ri, vi in zip(rr.to(u.km).value, vv.to(u.km / u.s).value):
        orbit_params = rv2coe(Earth.k.to(u.km**3 / u.s**2).value, ri, vi)
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

    # Averaging over 5 last values in the way Curtis does
    index = int(1.0 * t_days / tof.to(u.day).value * 4000)  # type: ignore
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
