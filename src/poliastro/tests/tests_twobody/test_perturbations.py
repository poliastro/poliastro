import pytest
from poliastro.twobody.propagation import cowell
import numpy as np
from poliastro.twobody.rv import rv2coe
from astropy import units as u
from poliastro.util import norm
from poliastro.twobody.perturbations import J2_perturbation, atmospheric_drag
from poliastro.bodies import Earth
from astropy.tests.helper import assert_quantity_allclose
from poliastro.twobody import Orbit


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
    # given the expression for \dot{a} / a, aproximate \Delta a \approx F_a * \Delta t

    R = Earth.R.to(u.km).value
    k = Earth.k.to(u.km**3 / u.s**2).value

    # parameters of a circular orbit with h = 10 km
    h = 250  # km
    r0 = np.array([R + h, 0, 0])  # km
    v0 = np.array([0, np.sqrt(k / (R + h)), 0])  # km/s

    # parameters of a body
    C_D = 2.2  # dimentionless
    A = ((np.pi / 4.0) * (u.m**2)).to(u.km**2).value  # km^2
    m = 100  # kg
    B = C_D * A / m

    # parameters of the atmosphere
    rho0 = Earth.rho0.to(u.kg / u.km**3).value  # kg/km^3
    H0 = Earth.H0.to(u.km).value

    orbit = Orbit.from_vectors(Earth, r0 * u.km, v0 * u.km / u.s)
    tof = 100000  # s

    dr_expected = -B * tof * rho0 * np.exp(-(norm(r0) - R) / H0) * np.sqrt(k * norm(r0))

    r, v = cowell(orbit, tof, ad=atmospheric_drag, R=R, B=B, H0=H0, rho0=rho0)

    assert_quantity_allclose(norm(r) - norm(r0), dr_expected, rtol=1e-2)
