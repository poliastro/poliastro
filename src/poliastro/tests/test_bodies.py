import math
import os

import pytest
import requests
import spiceypy as spice
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time

from poliastro import bodies, constants
from poliastro.bodies import (
    Earth,
    Jupiter,
    Mars,
    Mercury,
    Moon,
    Neptune,
    Pluto,
    Saturn,
    Sun,
    Uranus,
    Venus,
)


def test_body_has_k_given_in_constructor():
    k = 3.98e5 * u.km ** 3 / u.s ** 2
    earth = bodies.Body(None, k, "")
    assert earth.k == k


def test_body_from_parameters_raises_valueerror_if_k_units_not_correct():
    wrong_k = 4902.8 * u.kg
    _name = _symbol = ""
    _R = 0
    with pytest.raises(u.UnitsError) as excinfo:
        bodies.Body.from_parameters(None, wrong_k, _name, _symbol, _R)
    assert (
        "UnitsError: Argument 'k' to function 'from_parameters' must be in units convertible to 'km3 / s2'."
        in excinfo.exconly()
    )


def test_body_printing_has_name_and_symbol():
    name = "2 Pallas"
    symbol = u"\u26b4"
    k = 1.41e10 * u.m ** 3 / u.s ** 2
    pallas2 = bodies.Body(None, k, name, symbol)
    assert name in str(pallas2)
    assert symbol in str(pallas2)


def test_earth_has_k_given_in_literature():
    expected_k = 3.986004418e14 * u.m ** 3 / u.s ** 2
    k = bodies.Earth.k
    assert_quantity_allclose(k.decompose([u.km, u.s]), expected_k)


def test_body_kwargs():
    name = "2 Pallas"
    symbol = u"\u26b4"
    k = 1.41e10 * u.m ** 3 / u.s ** 2
    pallas2 = bodies.Body(None, k, name, symbol)
    assert pallas2.kwargs == {}
    pallas2 = bodies.Body(None, k, name, symbol, custom="data")
    assert "custom" in pallas2.kwargs


def test_from_relative():
    TRAPPIST1 = bodies.Body.from_relative(
        reference=bodies.Sun,
        parent=None,
        k=0.08,  # Relative to the Sun
        name="TRAPPIST",
        symbol=None,
        R=0.114,
    )  # Relative to the Sun

    # check values properly calculated
    VALUECHECK = bodies.Body.from_relative(
        reference=bodies.Earth,
        parent=TRAPPIST1,
        k=1,
        name="VALUECHECK",
        symbol=None,
        R=1,
    )
    assert bodies.Earth.k == VALUECHECK.k
    assert bodies.Earth.R == VALUECHECK.R


class TestRotElements:
    """Orientation models to calculate the rotational elements (ICRF right ascension
    and declination and prime meridian location) are given in the kernel file
    "pck00010.tpc", provided by Nasa. These models have been used to calculate
    the rotational elements.
    The kernel file provides the coefficients of polynomials and trigonometric terms of
    the models. Hence, in the tests, coefficients are taken from the kernel file and
    plugged into the polynomial/trigonometric terms.
    Url for kernel file is "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc"
    """

    # Download the kernel file
    kernel_folder = os.path.join(os.getcwd(), "src/poliastro/tests/kernels")
    kernel_name = os.path.join(kernel_folder, "pck00010.tpc")
    kernel_url = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc"
    r = requests.get(kernel_url)
    with open(kernel_name, "wb") as f:
        f.write(r.content)
    spice.furnsh(kernel_name)

    @pytest.fixture()
    def set_epoch(self):
        time = ["2000-01-01T00:00:00.00"]
        epoch = Time(time, format="isot", scale="utc")
        T = (epoch.tdb - constants.J2000).to("day").value / 36525
        d = (epoch.tdb - constants.J2000).to("day").value
        return epoch, T, d

    def test_Sun(self, set_epoch):
        epoch, T, d = set_epoch
        value_from_poliastro = [i.value for i in Sun.rot_elements_at_epoch(epoch)]

        body_code = spice.bodn2c("Sun")

        BODY_POLE_RA = spice.gdpool("BODY" + str(body_code) + "_POLE_RA", 0, 3)
        BODY_POLE_DEC = spice.gdpool("BODY" + str(body_code) + "_POLE_DEC", 0, 3)
        BODY_PM = spice.gdpool("BODY" + str(body_code) + "_PM", 0, 3)

        body_RA = BODY_POLE_RA[0]
        body_DEC = BODY_POLE_DEC[0]
        body_PM = BODY_PM[0] + BODY_PM[1] * d
        value_from_spice = [body_RA, body_DEC, body_PM]

        for i, j in zip(value_from_poliastro, value_from_spice):
            assert_quantity_allclose(i, j)

    def test_Mercury(self, set_epoch):
        epoch, T, d = set_epoch
        value_from_poliastro = [i.value for i in Mercury.rot_elements_at_epoch(epoch)]

        body_code = spice.bodn2c("Mercury")

        BODY_POLE_RA = spice.gdpool("BODY" + str(body_code) + "_POLE_RA", 0, 3)
        BODY_POLE_DEC = spice.gdpool("BODY" + str(body_code) + "_POLE_DEC", 0, 3)
        BODY_PM = spice.gdpool("BODY" + str(body_code) + "_PM", 0, 3)

        body_RA = BODY_POLE_RA[0] + BODY_POLE_RA[1] * T
        body_DEC = BODY_POLE_DEC[0] + BODY_POLE_DEC[1] * T
        body_PM = BODY_PM[0] + BODY_PM[1] * d
        value_from_spice = [body_RA, body_DEC, body_PM]

        for i, j in zip(value_from_poliastro, value_from_spice):
            assert_quantity_allclose(i, j)

    def test_Venus(self, set_epoch):
        epoch, T, d = set_epoch
        value_from_poliastro = [i.value for i in Venus.rot_elements_at_epoch(epoch)]

        body_code = spice.bodn2c("Venus")

        BODY_POLE_RA = spice.gdpool("BODY" + str(body_code) + "_POLE_RA", 0, 3)
        BODY_POLE_DEC = spice.gdpool("BODY" + str(body_code) + "_POLE_DEC", 0, 3)
        BODY_PM = spice.gdpool("BODY" + str(body_code) + "_PM", 0, 3)

        body_RA = BODY_POLE_RA[0] + BODY_POLE_RA[1] * T
        body_DEC = BODY_POLE_DEC[0] + BODY_POLE_DEC[1] * T
        body_PM = BODY_PM[0] + BODY_PM[1] * d
        value_from_spice = [body_RA, body_DEC, body_PM]

        for i, j in zip(value_from_poliastro, value_from_spice):
            assert_quantity_allclose(i, j)

    def test_Earth(self, set_epoch):
        epoch, T, d = set_epoch
        value_from_poliastro = [i.value for i in Earth.rot_elements_at_epoch(epoch)]

        body_code = spice.bodn2c("Earth")

        BODY_POLE_RA = spice.gdpool("BODY" + str(body_code) + "_POLE_RA", 0, 3)
        BODY_POLE_DEC = spice.gdpool("BODY" + str(body_code) + "_POLE_DEC", 0, 3)
        BODY_PM = spice.gdpool("BODY" + str(body_code) + "_PM", 0, 3)

        body_RA = BODY_POLE_RA[0] + BODY_POLE_RA[1] * T
        body_DEC = BODY_POLE_DEC[0] + BODY_POLE_DEC[1] * T
        body_PM = BODY_PM[0] + BODY_PM[1] * d
        value_from_spice = [body_RA, body_DEC, body_PM]

        for i, j in zip(value_from_poliastro, value_from_spice):
            assert_quantity_allclose(i, j)

    def test_Mars(self, set_epoch):
        epoch, T, d = set_epoch
        value_from_poliastro = [i.value for i in Mars.rot_elements_at_epoch(epoch)]

        body_code = spice.bodn2c("Mars")

        BODY_POLE_RA = spice.gdpool("BODY" + str(body_code) + "_POLE_RA", 0, 3)
        BODY_POLE_DEC = spice.gdpool("BODY" + str(body_code) + "_POLE_DEC", 0, 3)
        BODY_PM = spice.gdpool("BODY" + str(body_code) + "_PM", 0, 3)

        body_RA = BODY_POLE_RA[0] + BODY_POLE_RA[1] * T
        body_DEC = BODY_POLE_DEC[0] + BODY_POLE_DEC[1] * T
        body_PM = BODY_PM[0] + BODY_PM[1] * d
        value_from_spice = [body_RA, body_DEC, body_PM]

        for i, j in zip(value_from_poliastro, value_from_spice):
            assert_quantity_allclose(i, j)

    def test_Jupiter(self, set_epoch):
        epoch, T, d = set_epoch
        value_from_poliastro = [i.value for i in Jupiter.rot_elements_at_epoch(epoch)]

        body_code = spice.bodn2c("Jupiter")

        BODY_POLE_RA = spice.gdpool("BODY" + str(body_code) + "_POLE_RA", 0, 3)
        BODY_POLE_DEC = spice.gdpool("BODY" + str(body_code) + "_POLE_DEC", 0, 3)
        BODY_PM = spice.gdpool("BODY" + str(body_code) + "_PM", 0, 3)
        BODY5_NUT_PREC_ANGLES = spice.gdpool("BODY5_NUT_PREC_ANGLES", 0, 30)
        BODY599_NUT_PREC_RA = spice.gdpool("BODY599_NUT_PREC_RA", 0, 15)
        BODY599_NUT_PREC_DEC = spice.gdpool("BODY599_NUT_PREC_DEC", 0, 15)

        Ja = BODY5_NUT_PREC_ANGLES[20] + BODY5_NUT_PREC_ANGLES[21] * T
        Jb = BODY5_NUT_PREC_ANGLES[22] + BODY5_NUT_PREC_ANGLES[23] * T
        Jc = BODY5_NUT_PREC_ANGLES[24] + BODY5_NUT_PREC_ANGLES[25] * T
        Jd = BODY5_NUT_PREC_ANGLES[26] + BODY5_NUT_PREC_ANGLES[28] * T
        Je = BODY5_NUT_PREC_ANGLES[28] + BODY5_NUT_PREC_ANGLES[29] * T

        body_RA = (
            BODY_POLE_RA[0]
            + BODY_POLE_RA[1] * T
            + BODY599_NUT_PREC_RA[10] * math.sin(math.radians(Ja))
            + BODY599_NUT_PREC_RA[11] * math.sin(math.radians(Jb))
            + BODY599_NUT_PREC_RA[12] * math.sin(math.radians(Jc))
            + BODY599_NUT_PREC_RA[13] * math.sin(math.radians(Jd))
            + BODY599_NUT_PREC_RA[14] * math.sin(math.radians(Je))
        )

        body_DEC = (
            BODY_POLE_DEC[0]
            + BODY_POLE_DEC[1] * T
            + BODY599_NUT_PREC_DEC[10] * math.cos(math.radians(Ja))
            + BODY599_NUT_PREC_DEC[11] * math.cos(math.radians(Jb))
            + BODY599_NUT_PREC_DEC[12] * math.cos(math.radians(Jc))
            + BODY599_NUT_PREC_DEC[13] * math.cos(math.radians(Jd))
            + BODY599_NUT_PREC_DEC[14] * math.cos(math.radians(Je))
        )

        body_PM = BODY_PM[0] + BODY_PM[1] * d

        value_from_spice = [body_RA, body_DEC, body_PM]

        for i, j in zip(value_from_poliastro, value_from_spice):
            assert_quantity_allclose(i, j)

    def test_Saturn(self, set_epoch):
        epoch, T, d = set_epoch
        value_from_poliastro = [i.value for i in Saturn.rot_elements_at_epoch(epoch)]

        body_code = spice.bodn2c("Saturn")

        BODY_POLE_RA = spice.gdpool("BODY" + str(body_code) + "_POLE_RA", 0, 3)
        BODY_POLE_DEC = spice.gdpool("BODY" + str(body_code) + "_POLE_DEC", 0, 3)
        BODY_PM = spice.gdpool("BODY" + str(body_code) + "_PM", 0, 3)

        body_RA = BODY_POLE_RA[0] + BODY_POLE_RA[1] * T
        body_DEC = BODY_POLE_DEC[0] + BODY_POLE_DEC[1] * T
        body_PM = BODY_PM[0] + BODY_PM[1] * d
        value_from_spice = [body_RA, body_DEC, body_PM]

        for i, j in zip(value_from_poliastro, value_from_spice):
            assert_quantity_allclose(i, j)

    def test_Uranus(self, set_epoch):
        epoch, T, d = set_epoch
        value_from_poliastro = [i.value for i in Uranus.rot_elements_at_epoch(epoch)]

        body_code = spice.bodn2c("Uranus")

        BODY_POLE_RA = spice.gdpool("BODY" + str(body_code) + "_POLE_RA", 0, 3)
        BODY_POLE_DEC = spice.gdpool("BODY" + str(body_code) + "_POLE_DEC", 0, 3)
        BODY_PM = spice.gdpool("BODY" + str(body_code) + "_PM", 0, 3)

        body_RA = BODY_POLE_RA[0] + BODY_POLE_RA[1] * T
        body_DEC = BODY_POLE_DEC[0] + BODY_POLE_DEC[1] * T
        body_PM = BODY_PM[0] + BODY_PM[1] * d
        value_from_spice = [body_RA, body_DEC, body_PM]

        for i, j in zip(value_from_poliastro, value_from_spice):
            assert_quantity_allclose(i, j)

    def test_Neptune(self, set_epoch):
        epoch, T, d = set_epoch
        value_from_poliastro = [i.value for i in Neptune.rot_elements_at_epoch(epoch)]

        body_code = spice.bodn2c("Neptune")

        BODY_POLE_RA = spice.gdpool("BODY" + str(body_code) + "_POLE_RA", 0, 3)
        BODY_POLE_DEC = spice.gdpool("BODY" + str(body_code) + "_POLE_DEC", 0, 3)
        BODY_PM = spice.gdpool("BODY" + str(body_code) + "_PM", 0, 3)

        BODY899_NUT_PREC_RA = spice.gdpool("BODY899_NUT_PREC_RA", 0, 8)
        BODY899_NUT_PREC_DEC = spice.gdpool("BODY899_NUT_PREC_DEC", 0, 8)
        BODY899_NUT_PREC_PM = spice.gdpool("BODY899_NUT_PREC_PM", 0, 8)
        BODY8_NUT_PREC_ANGLES = spice.gdpool("BODY8_NUT_PREC_ANGLES", 0, 35)

        N = BODY8_NUT_PREC_ANGLES[0] + BODY8_NUT_PREC_ANGLES[1] * T
        body_RA = (
            BODY_POLE_RA[0]
            + BODY_POLE_RA[1] * T
            + BODY899_NUT_PREC_RA[0] * math.sin(math.radians(N))
        )
        body_DEC = (
            BODY_POLE_DEC[0]
            + BODY_POLE_DEC[1] * T
            + BODY899_NUT_PREC_DEC[0] * math.cos(math.radians(N))
        )
        body_PM = (
            BODY_PM[0]
            + BODY_PM[1] * d
            + BODY899_NUT_PREC_PM[0] * math.sin(math.radians(N))
        )
        value_from_spice = [body_RA, body_DEC, body_PM]

        for i, j in zip(value_from_poliastro, value_from_spice):
            assert_quantity_allclose(i, j)

    def test_Pluto(self, set_epoch):
        epoch, T, d = set_epoch
        value_from_poliastro = [i.value for i in Pluto.rot_elements_at_epoch(epoch)]

        body_code = spice.bodn2c("Pluto")

        BODY_POLE_RA = spice.gdpool("BODY" + str(body_code) + "_POLE_RA", 0, 3)
        BODY_POLE_DEC = spice.gdpool("BODY" + str(body_code) + "_POLE_DEC", 0, 3)
        BODY_PM = spice.gdpool("BODY" + str(body_code) + "_PM", 0, 3)

        body_RA = BODY_POLE_RA[0] + BODY_POLE_RA[1] * T
        body_DEC = BODY_POLE_DEC[0] + BODY_POLE_DEC[1] * T
        body_PM = BODY_PM[0] + BODY_PM[1] * d
        value_from_spice = [body_RA, body_DEC, body_PM]

        for i, j in zip(value_from_poliastro, value_from_spice):
            assert_quantity_allclose(i, j)

    def test_Moon(self, set_epoch):
        epoch, T, d = set_epoch
        value_from_poliastro = [i.value for i in Moon.rot_elements_at_epoch(epoch)]

        body_code = spice.bodn2c("Moon")

        BODY_POLE_RA = spice.gdpool("BODY" + str(body_code) + "_POLE_RA", 0, 3)
        BODY_POLE_DEC = spice.gdpool("BODY" + str(body_code) + "_POLE_DEC", 0, 3)
        BODY_PM = spice.gdpool("BODY" + str(body_code) + "_PM", 0, 3)
        BODY3_NUT_PREC_ANGLES = spice.gdpool("BODY3_NUT_PREC_ANGLES", 0, 30)
        BODY301_NUT_PREC_RA = spice.gdpool("BODY301_NUT_PREC_RA", 0, 15)
        BODY301_NUT_PREC_DEC = spice.gdpool("BODY301_NUT_PREC_DEC", 0, 15)
        BODY301_NUT_PREC_PM = spice.gdpool("BODY301_NUT_PREC_PM", 0, 15)

        E1 = BODY3_NUT_PREC_ANGLES[0] + BODY3_NUT_PREC_ANGLES[1] * T
        E2 = BODY3_NUT_PREC_ANGLES[2] + BODY3_NUT_PREC_ANGLES[3] * T
        E3 = BODY3_NUT_PREC_ANGLES[4] + BODY3_NUT_PREC_ANGLES[5] * T
        E4 = BODY3_NUT_PREC_ANGLES[6] + BODY3_NUT_PREC_ANGLES[7] * T
        E5 = BODY3_NUT_PREC_ANGLES[8] + BODY3_NUT_PREC_ANGLES[9] * T
        E6 = BODY3_NUT_PREC_ANGLES[10] + BODY3_NUT_PREC_ANGLES[11] * T
        E7 = BODY3_NUT_PREC_ANGLES[12] + BODY3_NUT_PREC_ANGLES[13] * T
        E8 = BODY3_NUT_PREC_ANGLES[14] + BODY3_NUT_PREC_ANGLES[15] * T
        E9 = BODY3_NUT_PREC_ANGLES[16] + BODY3_NUT_PREC_ANGLES[17] * T
        E10 = BODY3_NUT_PREC_ANGLES[18] + BODY3_NUT_PREC_ANGLES[19] * T
        E11 = BODY3_NUT_PREC_ANGLES[20] + BODY3_NUT_PREC_ANGLES[21] * T
        E12 = BODY3_NUT_PREC_ANGLES[22] + BODY3_NUT_PREC_ANGLES[23] * T
        E13 = BODY3_NUT_PREC_ANGLES[24] + BODY3_NUT_PREC_ANGLES[25] * T
        body_RA = (
            BODY_POLE_RA[0]
            + BODY_POLE_RA[1] * T
            + BODY301_NUT_PREC_RA[0] * math.sin(math.radians(E1))
            + BODY301_NUT_PREC_RA[1] * math.sin(math.radians(E2))
            + BODY301_NUT_PREC_RA[2] * math.sin(math.radians(E3))
            + BODY301_NUT_PREC_RA[3] * math.sin(math.radians(E4))
            + BODY301_NUT_PREC_RA[4] * math.sin(math.radians(E5))
            + BODY301_NUT_PREC_RA[5] * math.sin(math.radians(E6))
            + BODY301_NUT_PREC_RA[6] * math.sin(math.radians(E7))
            + BODY301_NUT_PREC_RA[7] * math.sin(math.radians(E8))
            + BODY301_NUT_PREC_RA[8] * math.sin(math.radians(E9))
            + BODY301_NUT_PREC_RA[9] * math.sin(math.radians(E10))
            + BODY301_NUT_PREC_RA[10] * math.sin(math.radians(E11))
            + BODY301_NUT_PREC_RA[11] * math.sin(math.radians(E12))
            + BODY301_NUT_PREC_RA[12] * math.sin(math.radians(E13))
        )
        body_DEC = (
            BODY_POLE_DEC[0]
            + BODY_POLE_DEC[1] * T
            + BODY301_NUT_PREC_DEC[0] * math.cos(math.radians(E1))
            + BODY301_NUT_PREC_DEC[1] * math.cos(math.radians(E2))
            + BODY301_NUT_PREC_DEC[2] * math.cos(math.radians(E3))
            + BODY301_NUT_PREC_DEC[3] * math.cos(math.radians(E4))
            + BODY301_NUT_PREC_DEC[4] * math.cos(math.radians(E5))
            + BODY301_NUT_PREC_DEC[5] * math.cos(math.radians(E6))
            + BODY301_NUT_PREC_DEC[6] * math.cos(math.radians(E7))
            + BODY301_NUT_PREC_DEC[7] * math.cos(math.radians(E8))
            + BODY301_NUT_PREC_DEC[8] * math.cos(math.radians(E9))
            + BODY301_NUT_PREC_DEC[9] * math.cos(math.radians(E10))
            + BODY301_NUT_PREC_DEC[10] * math.cos(math.radians(E11))
            + BODY301_NUT_PREC_DEC[11] * math.cos(math.radians(E12))
            + BODY301_NUT_PREC_DEC[12] * math.cos(math.radians(E13))
        )
        body_PM = (
            BODY_PM[0]
            + BODY_PM[1] * d
            + BODY_PM[2] * d ** 2
            + BODY301_NUT_PREC_PM[0] * math.sin(math.radians(E1))
            + BODY301_NUT_PREC_PM[1] * math.sin(math.radians(E2))
            + BODY301_NUT_PREC_PM[2] * math.sin(math.radians(E3))
            + BODY301_NUT_PREC_PM[3] * math.sin(math.radians(E4))
            + BODY301_NUT_PREC_PM[4] * math.sin(math.radians(E5))
            + BODY301_NUT_PREC_PM[5] * math.sin(math.radians(E6))
            + BODY301_NUT_PREC_PM[6] * math.sin(math.radians(E7))
            + BODY301_NUT_PREC_PM[7] * math.sin(math.radians(E8))
            + BODY301_NUT_PREC_PM[8] * math.sin(math.radians(E9))
            + BODY301_NUT_PREC_PM[9] * math.sin(math.radians(E10))
            + BODY301_NUT_PREC_PM[10] * math.sin(math.radians(E11))
            + BODY301_NUT_PREC_PM[11] * math.sin(math.radians(E12))
            + BODY301_NUT_PREC_PM[12] * math.sin(math.radians(E13))
        )
        value_from_spice = [body_RA, body_DEC, body_PM]

        for i, j in zip(value_from_poliastro, value_from_spice):
            assert_quantity_allclose(i, j)

    os.remove(kernel_name)
