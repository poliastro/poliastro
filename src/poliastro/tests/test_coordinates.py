import astropy.units as u
from astropy import time
from astropy.tests.helper import assert_quantity_allclose
from numpy.linalg import norm

from poliastro import bodies, coordinates
from poliastro.examples import molniya
from poliastro.twobody.orbit import Orbit


# Note that function are tested using astropy current builtin ephemeris.
# Horizons uses JPL ephemeris DE431, so expected values are hardcoded,
# instead of being obtained using Horizons.
def test_body_centered_to_icrs_transformation():

    vexpress_r_venus = [
        -2.707041558060933e03,
        1.112962479175306e04,
        -3.792944408664889e04,
    ] * u.km
    vexpress_v_venus = (
        [-2.045118200275925e-01, 7.978578482960554e-01, 2.664944903217139e00]
        * u.km
        / u.s
    )

    expected_r = [
        -3.47202219448080286e07,
        9.16853879708216339e07,
        4.34117810525591150e07,
    ] * u.km
    expected_v = (
        [-3.34053728321152121e01, -1.16604776013667291e01, -2.39943678872506838e-01]
        * u.km
        / u.s
    )

    r, v = coordinates.body_centered_to_icrs(
        vexpress_r_venus,
        vexpress_v_venus,
        bodies.Venus,
        time.Time("2014-08-23 00:00", scale="tdb"),
    )

    assert_quantity_allclose(r, expected_r)
    assert_quantity_allclose(v, expected_v)


def test_icrs_to_body_centered_transformation():
    vexpress_r_icrs = [
        -3.472125578094885e07,
        9.168528034176737e07,
        4.341160627674723e07,
    ] * u.km
    vexpress_v_icrs = (
        [-3.340574196483147e01, -1.165974037637970e01, -2.395829145441408e-01]
        * u.km
        / u.s
    )

    expected_r = [
        -3.74486105008138566e03,
        1.10085874027602295e04,
        -3.80681106516677464e04,
    ] * u.km
    expected_v = (
        [-2.04845025352488774e-01, 7.98692896032012989e-01, 2.66498465286454023e00]
        * u.km
        / u.s
    )

    r, v = coordinates.icrs_to_body_centered(
        vexpress_r_icrs,
        vexpress_v_icrs,
        bodies.Venus,
        time.Time("2014-08-23 00:00", scale="tdb"),
    )

    assert_quantity_allclose(r, expected_r)
    assert_quantity_allclose(v, expected_v)


def test_inertial_body_centered_to_pqw():
    molniya_r_peri, molniya_v_peri = coordinates.inertial_body_centered_to_pqw(
        molniya.r, molniya.v, bodies.Earth
    )

    molniya_peri = Orbit.from_vectors(
        bodies.Earth, molniya_r_peri, molniya_v_peri, molniya.epoch
    )

    assert_quantity_allclose(molniya_peri.e_vec[-2:], [0, 0], atol=1e-12)
    assert_quantity_allclose(norm(molniya_peri.e_vec), norm(molniya.e_vec))
