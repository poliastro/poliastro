import numpy as np
import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from poliastro.atmosphere.models import COESA62, COESA76
from poliastro.atmosphere.util import (
    geometric_to_geopotential,
    geopotential_to_geometric,
)
from poliastro.bodies import Earth

# U.S Standard Atmosphere base layer parameters

# Lower layers
Hb_lower62 = [0, 11, 20, 32, 47, 52, 61, 79] * u.km
Zb_lower62 = [
    geopotential_to_geometric(H, Earth.R).to(u.km).value for H in Hb_lower62
] * u.km

# Upper layers
Zb_upper62 = [
    90,
    100,
    110,
    120,
    150,
    160,
    170,
    190,
    230,
    300,
    400,
    500,
    600,
    700,
] * u.km
Hb_upper62 = [
    geometric_to_geopotential(Z, Earth.R).to(u.km).value for Z in Zb_upper62
] * u.km

# Layers
Hb_table62 = np.append(Hb_lower62, Hb_upper62).value * u.km
Zb_table62 = np.append(Zb_lower62, Zb_upper62).value * u.km
Tb_table62 = [
    288.15,
    216.65,
    216.65,
    228.65,
    270.65,
    270.65,
    252.65,
    180.65,
    180.65,
    210.65,
    260.65,
    360.65,
    960.65,
    1110.65,
    1210.65,
    1350.65,
    1550.65,
    1830.65,
    2160.65,
    2420.65,
    2590.65,
    2700.65,
] * u.K
Lb_table62 = (
    [
        -6.5,
        0.0,
        1.0,
        2.8,
        0.0,
        -2.0,
        -4.0,
        0.0,
        3.0,
        5.0,
        10.0,
        20.0,
        15.0,
        10.0,
        7.0,
        5.0,
        4.0,
        3.3,
        2.6,
        1.7,
        1.1,
        0.0,
    ]
    * u.K
    / u.km
)

# Pressure table
pb_table62 = [
    1.01325e3,
    2.2632e2,
    5.47487e1,
    8.68014e0,
    1.10905e0,
    5.90005e-1,
    1.82099e-1,
    1.0377e-2,
    1.6438e-3,
    3.0075e-4,
    7.3544e-5,
    2.5217e-5,
    5.0617e-6,
    3.6943e-6,
    2.7926e-6,
    1.6852e-6,
    6.9604e-7,
    1.8838e-7,
    4.0304e-8,
    1.0957e-8,
    3.4502e-9,
    1.1918e-9,
] * u.mbar


# U.S Standard Atmosphere 1976
# Lower layers
Hb_lower76 = [0, 11, 20, 32, 47, 51, 71] * u.km
Zb_lower76 = [
    geopotential_to_geometric(H, Earth.R).to(u.km).value for H in Hb_lower76
] * u.km

# Upper layers
Zb_upper76 = [86, 91, 110, 120, 500, 1000] * u.km
Hb_upper76 = [
    geometric_to_geopotential(Z, Earth.R).to(u.km).value for Z in Zb_upper76
] * u.km

# Geopotential and geometric layers
Hb_table76 = np.append(Hb_lower76, Hb_upper76).value * u.km
Zb_table76 = np.append(Zb_lower76, Zb_upper76).value * u.km

# Temperature table
Tb_table76 = [
    288.15,
    216.65,
    216.650,
    228.65,
    270.65,
    270.65,
    214.65,
    186.95,
    187.36,
    254.93,
    397.91,
    2019.69,
    7351.15,
] * u.K

# Gradient temperature
Lb_table76 = (
    [-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0, 0.0, 12.0, 0.0, 0.0, 0.0] * u.K / u.km
)

# Pressure table
pb_table76 = [
    1.01325e3,
    2.2632e2,
    5.4748e1,
    8.6801e0,
    1.1090e0,
    6.6938e-1,
    3.9564e-2,
    3.7338e-3,
    1.5381e-3,
    7.1042e-5,
    2.5382e-5,
    3.0236e-9,
    7.5138e-11,
] * u.mbar


# Main atmosphere instances
table_coesa62 = [Zb_table62, Hb_table62, Tb_table62, Lb_table62, pb_table62]
table_coesa76 = [Zb_table76, Hb_table76, Tb_table76, Lb_table76, pb_table76]
coesa62 = COESA62()
coesa76 = COESA76()

# Pytest parametrization
atm_list = [coesa62, coesa76]
atm_dic = {coesa62: table_coesa62, coesa76: table_coesa76}


@pytest.mark.parametrize("atm", atm_list)
def test_base_layer_geometric(atm):
    for Z, Zb in zip(atm_dic[atm][0], atm._Zb_table):
        assert_quantity_allclose(Zb, Z)


@pytest.mark.parametrize("atm", atm_list)
def test_base_layer_geopotential(atm):
    for H, Hb in zip(atm_dic[atm][1], atm._Hb_table):
        assert_quantity_allclose(Hb, H)


@pytest.mark.parametrize("atm", atm_list)
def test_base_layer_temperature(atm):
    for T, Tb in zip(atm_dic[atm][2], atm._Tb_table):
        assert_quantity_allclose(Tb, T)


@pytest.mark.parametrize("atm", atm_list)
def test_base_layer_gradient(atm):
    for L, Lb in zip(atm_dic[atm][3], atm._Lb_table):
        assert_quantity_allclose(Lb, L)


@pytest.mark.parametrize("atm", atm_list)
def test_pressure_tables(atm):
    for p, pb in zip(atm_dic[atm][4], atm._pb_table):
        assert_quantity_allclose(pb, p)


@pytest.mark.parametrize("table", table_coesa62)
def test_base_layer_length_coesa62(table):
    assert len(table) == 22


@pytest.mark.parametrize("table", table_coesa76)
def test_base_layer_length_coesa76(table):
    assert len(table) == 13


def test_not_implemented_pressure_coesa62():
    with pytest.raises(NotImplementedError) as excinfo:
        coesa62.pressure(100 * u.km)
    assert (
        "Pressure in COESA62 has just been implemented up to 90km." in excinfo.exconly()
    )


def test_not_implemented_pressure_coesa76():
    with pytest.raises(NotImplementedError) as excinfo:
        coesa76.pressure(100 * u.km)
    assert (
        "Pressure in COESA76 has just been implemented up to 86km." in excinfo.exconly()
    )


def test_not_implemented_density_coesa62():
    with pytest.raises(NotImplementedError) as excinfo:
        coesa62.density(100 * u.km)
    assert (
        "Density in COESA62 has just been implemented up to 90km." in excinfo.exconly()
    )


def test_not_implemented_density_coesa76():
    with pytest.raises(NotImplementedError) as excinfo:
        coesa76.density(100 * u.km)
    assert (
        "Density in COESA76 has just been implemented up to 86km." in excinfo.exconly()
    )


def test_limit_altitude_coesa62():
    with pytest.raises(ValueError) as excinfo:
        coesa62.temperature(901 * u.km)
    assert (
        "Geometric altitude must be in range [0.0 km, 700.0 km]." in excinfo.exconly()
    )


def test_limit_altitude_coesa76():
    with pytest.raises(ValueError) as excinfo:
        coesa76.temperature(1001 * u.km)
    assert (
        "Geometric altitude must be in range [0.0 km, 1000.0 km]." in excinfo.exconly()
    )
