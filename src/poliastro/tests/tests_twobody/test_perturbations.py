import pytest
import functools
import numpy as np

from scipy.integrate import solve_ivp
from poliastro.integrators import DOP835

from astropy.time import Time
from poliastro.twobody.propagation import cowell, func_twobody
from poliastro.twobody.rv import rv2coe
from astropy import units as u
from poliastro.util import norm
from poliastro.twobody.perturbations import J2_perturbation, atmospheric_drag, third_body
from poliastro.bodies import Earth, Moon
from astropy.tests.helper import assert_quantity_allclose
from poliastro.twobody import Orbit
from poliastro.coordinates import transform
from astropy.coordinates import ICRS, GCRS


def test_J2_propagation_Earth():
    # from Curtis example 12.2:
    r0 = np.array([-2384.46, 5729.01, 3050.46])  # km
    v0 = np.array([-7.36138, -2.98997, 1.64354])  # km/s
    k = Earth.k.to(u.km**3 / u.s**2).value

    orbit = Orbit.from_vectors(Earth, r0 * u.km, v0 * u.km / u.s)

    tof = (48.0 * u.h).to(u.s).value
    r, v = cowell(orbit, tof, ad=J2_perturbation, J2=Earth.J2.value, R=Earth.R.to(u.km).value)

    _, _, _, raan0, argp0, _ = rv2coe(k, r0, v0)
    _, _, _, raan, argp, _ = rv2coe(k, r, v)

    raan_variation_rate = (raan - raan0) / tof
    argp_variation_rate = (argp - argp0) / tof

    raan_variation_rate = (raan_variation_rate * u.rad / u.s).to(u.deg / u.h)
    argp_variation_rate = (argp_variation_rate * u.rad / u.s).to(u.deg / u.h)

    assert_quantity_allclose(raan_variation_rate, -0.172 * u.deg / u.h, rtol=1e-2)
    assert_quantity_allclose(argp_variation_rate, 0.282 * u.deg / u.h, rtol=1e-2)


def test_atmospheric_drag():
    # http://farside.ph.utexas.edu/teaching/celestial/Celestialhtml/node94.html#sair (10.148)
    # given the expression for \dot{r} / r, aproximate \Delta r \approx F_r * \Delta t

    R = Earth.R.to(u.km).value
    k = Earth.k.to(u.km**3 / u.s**2).value

    # parameters of a circular orbit with h = 250 km (any value would do, but not too small)
    orbit = Orbit.circular(Earth, 250 * u.km)
    r0, _ = orbit.rv()
    r0 = r0.to(u.km).value

    # parameters of a body
    C_D = 2.2  # dimentionless (any value would do)
    A = ((np.pi / 4.0) * (u.m**2)).to(u.km**2).value  # km^2
    m = 100  # kg
    B = C_D * A / m

    # parameters of the atmosphere
    rho0 = Earth.rho0.to(u.kg / u.km**3).value  # kg/km^3
    H0 = Earth.H0.to(u.km).value
    tof = 100000  # s

    dr_expected = -B * rho0 * np.exp(-(norm(r0) - R) / H0) * np.sqrt(k * norm(r0)) * tof
    # assuming the atmospheric decay during tof is small,
    # dr_expected = F_r * tof (Newton's integration formula), where
    # F_r = -B rho(r) |r|^2 sqrt(k / |r|^3) = -B rho(r) sqrt(k |r|)

    r, v = cowell(orbit, tof, ad=atmospheric_drag, R=R, C_D=C_D, A=A, m=m, H0=H0, rho0=rho0)

    assert_quantity_allclose(norm(r) - norm(r0), dr_expected, rtol=1e-2)


def test_cowell_works_with_small_perturbations():
    r0 = [-2384.46, 5729.01, 3050.46] * u.km
    v0 = [-7.36138, -2.98997, 1.64354] * u.km / u.s

    r_expected = [13179.39566663877121754922, -13026.25123408228319021873, -9852.66213692844394245185] * u.km
    v_expected = [2.78170542314378943516, 3.21596786944631274352, 0.16327165546278937791] * u.km / u.s

    initial = Orbit.from_vectors(Earth, r0, v0)

    def accel(t0, state, k):
        v_vec = state[3:]
        norm_v = (v_vec * v_vec).sum() ** .5
        return 1e-5 * v_vec / norm_v

    final = initial.propagate(3 * u.day, method=cowell, ad=accel)

    assert_quantity_allclose(final.r, r_expected)
    assert_quantity_allclose(final.v, v_expected)


def test_cowell_converges_with_small_perturbations():
    r0 = [-2384.46, 5729.01, 3050.46] * u.km
    v0 = [-7.36138, -2.98997, 1.64354] * u.km / u.s

    initial = Orbit.from_vectors(Earth, r0, v0)

    def accel(t0, state, k):
        v_vec = state[3:]
        norm_v = (v_vec * v_vec).sum() ** .5
        return 0.0 * v_vec / norm_v

    final = initial.propagate(initial.period, method=cowell, ad=accel)
    assert_quantity_allclose(final.r, initial.r)
    assert_quantity_allclose(final.v, initial.v)


def test_3rd_body():
    # example 12.11 from Howard Curtis
    epoch = Time(2454283.0, format='jd', scale='tdb')

    # propagate Moon in the interval (0, tof) and write it into the OdeSolution object
    tof = (60 * u.day).to(u.s).value

    moon = Orbit.from_body_ephem(Moon, epoch)
    moon = transform(moon, ICRS, GCRS)

    x, y, z = moon.r.to(u.km).value
    vx, vy, vz = moon.v.to(u.km / u.s).value

    u0 = np.array([x, y, z, vx, vy, vz])

    # Set the non Keplerian acceleration
    ad = lambda t0, u_, k_: (0, 0, 0)
    f_with_ad = functools.partial(func_twobody, k=Earth.k.to(u.km ** 3 / u.s ** 2).value, ad=ad, ad_kwargs=dict())

    moon_dense = solve_ivp(f_with_ad, (0, tof), u0,
                           rtol=1e-8, atol=1e-12, method=DOP835,
                           dense_output=True)

    ecc = 0.741 * u.one
    a = 26553.4 * u.km
    raan = 0.0 * u.deg
    inc = 63.4 * u.deg
    argp = 270 * u.deg
    nu = 0.0 * u.rad

    initial = Orbit.from_classical(Earth, a, ecc, inc, raan, argp, nu, epoch=epoch)

    r, v = cowell(initial, tof, rtol=1e-8, ad=third_body,
                  k_third=Moon.k.to(u.km**3 / u.s**2).value, third_body=moon_dense)
    _, _, inc_f, raan_f, argp_f, _ = rv2coe(Earth.k.to(u.km**3 / u.s**2).value, r, v)

    assert_quantity_allclose((raan_f * u.rad).to(u.deg) - 360 * u.deg, -0.06 * u.deg, rtol=1e-2)
    # assert_quantity_allclose((argp_f * u.rad - raan * u.rad).to(u.deg) - 360 * u.deg, -0.06 * u.deg, rtol=1e-2)
    # assert_quantity_allclose((raan_f * u.rad - raan * u.rad).to(u.deg) - 360 * u.deg, -0.06 * u.deg, rtol=1e-2)


def test_moon_at_right_position():
    # based on Cowell, example 12.10
    epoch = Time(2456498.8333, format='jd', scale='tdb')
    moon = Orbit.from_body_ephem(Moon, epoch)
    moon = transform(moon, ICRS, GCRS)
    r_expected = np.array([340883.595, -137975.741, -27888.668]) * u.km

    assert_quantity_allclose(moon.r, r_expected, rtol=1e-2)


def test_moon_propagates_right():
    # based on Cowell, example 12.10
    epoch = Time(2456498.8333, format='jd', scale='tdb')
    moon = Orbit.from_body_ephem(Moon, epoch)
    moon_ini = transform(moon, ICRS, GCRS)

    moon_final = moon_ini.propagate(20 * u.day)
    moon_check = transform(Orbit.from_body_ephem(Moon, Time(2456498.8333 + 20, format='jd', scale='tdb')), ICRS, GCRS)

    assert_quantity_allclose(moon_final.r, moon_check.r, rtol=1e-2)
