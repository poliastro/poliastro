""" Holds the Atmosphere class. """

import numpy as np
from astropy import units as u

from poliastro.atmosphere.util import (
    R_air,
    geometric_to_geopotential,
    geopotential_to_geometric,
    gravity,
)
from poliastro.constants import R_earth


class COESA(object):
    """ Class for U.S Standard Atmosphere models. """

    def __init__(self, table):
        """ Constructor for Atmosphere instances.

        Parameters
        ----------
        table: list 1x5
            List containing [Zb_table, Hb_table, Lb_table Tb_table, pb_table]
        """

        self._Zb_table = table[0]
        self._Hb_table = table[1]
        self._Tb_table = table[2]
        self._Lb_table = table[3]
        self._pb_table = table[4]

    def _get_heights(self, alt, geometric=True):
        """ Returns geometric and geopotential heights being passed altitude.

        Parameters
        ----------
        alt: ~astropy.units.Quantity
            Geometric/Geopotential height.
        geometric: boolean
            If `True`, assumes that `alt` argument is geometric kind.

        Returns
        -------
        Z: ~ astropy.units.Quantity
            Geometrical height.
        H: ~astropy.units.Quantity
            Geopotential height.
        """

        if geometric and (alt < self._Zb_table[0] or alt > self._Zb_table[-1]):
            raise ValueError(
                "Geometric altitude must be in range [{}, {}].".format(
                    self._Zb_table[0], self._Zb_table[-1]
                )
            )

        if not geometric and (alt < self._Hb_table[0] or alt > self._Hb_table[-1]):
            raise ValueError(
                "Geopotential altitude must be in range [{}, {}].".format(
                    self._Hb_table[0], self._Hb_table[-1]
                )
            )

        # Always work with geopotential height
        if geometric and alt in self._Zb_table:
            i = np.where(self._Zb_table == alt)
            Z = self._Zb_table[i]
            H = self._Hb_table[i]
        elif geometric:
            Z = alt
            H = geometric_to_geopotential(alt, R_earth)
        else:
            Z = geopotential_to_geometric(alt, R_earth)
            H = alt

        return Z, H

    def _get_index(self, a, n_set):
        """ Gets index inside array. """

        for i, b in enumerate(n_set):
            if i < len(n_set) and a > b:
                pass
            elif a == b:
                return i
            else:
                return i - 1

    def _get_base_index(self, H):
        """ Gets corresponding base coefficient for given geopotential altitude.

        Parameters
        ----------
        H: ~astropy.units.Quantity
            Geopotential altitude.

        Returns
        -------
        i: int
            Index position in table.
        """

        return self._get_index(H, self._Hb_table)

    def _get_base_parameters(self, H):
        """ Returns base layer parameters.

        Parameters
        ----------
        H: ~astropy.units.Quantity
            Geopotential height.

        Returns
        -------
        Zb: ~astropy.units.Quantity
            Geometric base height.
        Hb: ~astropy.units.Quantity
            Geopotential base height.
        Tb: ~astropy.units.Quantity
            Temperature at layer base.
        Lb: ~astropy.units.Quantity
            Temperature gradient at layer base.
        pb: ~lumped parameter model githubastropy.units.Quantity
            Pressure at layer base.
        """

        i = self._get_base_index(H)
        Zb = self._Zb_table[i]
        Hb = self._Hb_table[i]
        Tb = self._Tb_table[i]
        Lb = self._Lb_table[i]
        pb = self._pb_table[i]
        return Zb, Hb, Tb, Lb, pb


# COESA62 TABLES
# Lower layers (Z < 90km)
Hb_lower62 = [0.0, 11.0, 20.0, 32.0, 47.0, 52.0, 61.0, 79.0] * u.km
Zb_lower62 = [
    geopotential_to_geometric(H, R_earth).to(u.km).value for H in Hb_lower62
] * u.km

# Upper layers (Z => 90km)
Zb_upper62 = [
    90.0,
    100.0,
    110.0,
    120.0,
    150.0,
    160.0,
    170.0,
    190.0,
    230.0,
    300.0,
    400.0,
    500.0,
    600.0,
    700.0,
] * u.km
Hb_upper62 = [
    geometric_to_geopotential(Z, R_earth).to(u.km).value for Z in Zb_upper62
] * u.km

# Geopotential and geometric layers
Hb_table62 = np.append(Hb_lower62, Hb_upper62).value * u.km
Zb_table62 = np.append(Zb_lower62, Zb_upper62).value * u.km

# Temperature layers (defined by the molecular-scale temperature at each layer)
Tb_table62 = [
    288.15,
    216.650,
    216.650,
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

# Temperature gradient
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

# Pressure layers
pb_table62 = [
    1.01325e3,
    2.2632e2,
    5.47487e1,
    8.68014e0,
    1.109050,
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

# Generate table
TABLE_COESA62 = [Zb_table62, Hb_table62, Tb_table62, Lb_table62, pb_table62]


class COESA62(COESA):
    """ U.S Standard Atmosphere 1962.

    +--------+---------+---------+-----------+---------------+---------------+
    | Z (km) |  H (km) |  T (K)  |  p (mbar) | rho (kg / m3) | beta (K / km) |
    +--------+---------+---------+-----------+---------------+---------------+
    |   0.0  |   0.0   | 288.150 | 1.01325e3 |     1.2250    |      -6.5     |
    +--------+---------+---------+-----------+---------------+---------------+
    | 11.019 |   11.0  | 216.650 |  2.2632e2 |   3.6392e-1   |      0.0      |
    +--------+---------+---------+-----------+---------------+---------------+
    | 20.063 |   20.0  | 216.650 |  5.4749e1 |   8.8035e-2   |      1.0      |
    +--------+---------+---------+-----------+---------------+---------------+
    | 32.162 |   32.0  | 228.650 |  8.68014  |   1.3225e-2   |      2.8      |
    +--------+---------+---------+-----------+---------------+---------------+
    | 47.350 |   47.0  | 270.650 |  1.109050 |   1.4275e-3   |      0.0      |
    +--------+---------+---------+-----------+---------------+---------------+
    | 51.413 |   52.0  | 270.650 | 5.90005e-1|   7.5943e-4   |      -2.8     |
    +--------+---------+---------+-----------+---------------+---------------+
    | 61.591 |   61.0  | 252.650 | 1.82099e-1|   2.5109e-4   |      -2.0     |
    +--------+---------+---------+-----------+---------------+---------------+
    | 79.994 |   79.0  | 180.650 | 1.0377e-2 |    2.001e-5   |      0.0      |
    +--------+---------+---------+-----------+---------------+---------------+
    |  90.0  |  88.743 | 180.650 | 1.6438e-3 |    3.170e-6   |      0.0      |
    +--------+---------+---------+-----------+---------------+---------------+
    |  100.0 |  98.451 | 210.020 | 3.0075e-4 |    4.974e-7   |      5.0      |
    +--------+---------+---------+-----------+---------------+---------------+
    |  110.0 | 108.129 | 257.000 | 7.3544e-5 |    9.829e-8   |     10.0      |
    +--------+---------+---------+-----------+---------------+---------------+
    | 120.0  | 117.776 | 349.490 | 2.5217e-5 |    2.436e-8   |     20.0      |
    +--------+---------+---------+-----------+---------------+---------------+
    | 150.0  | 146.541 | 892.790 | 5.0617e-6 |    1.836e-9   |     15.0      |
    +--------+---------+---------+-----------+---------------+---------------+
    | 160.0  | 156.071 | 1022.23 | 3.6943e-6 |    1.159e-9   |     10.0      |
    +--------+---------+---------+-----------+---------------+---------------+
    | 170.0  | 165.571 | 1105.51 | 2.7926e-6 |    8.036e-10  |      7.0      |
    +--------+---------+---------+-----------+---------------+---------------+
    | 190.0  | 184.485 | 1205.50 | 1.6852e-6 |    4.347e-10  |      5.0      |
    +--------+---------+---------+-----------+---------------+---------------+
    | 230.0  | 221.967 | 1321.70 | 6.9604e-7 |    1.564e-10  |      4.0      |
    +--------+---------+---------+-----------+---------------+---------------+
    | 300.0  | 286.476 | 1432.11 | 1.8838e-7 |    3.585e-11  |      3.3      |
    +--------+---------+---------+-----------+---------------+---------------+
    | 400.0  | 376.312 | 1487.38 | 4.0304e-8 |    6.498e-12  |      2.6      |
    +--------+---------+---------+-----------+---------------+---------------+
    | 500.0  | 463.526 | 1499.22 | 1.0957e-8 |    1.577e-12  |      1.7      |
    +--------+---------+---------+-----------+---------------+---------------+
    | 600.0  | 548.230 | 1506.13 | 3.4502e-9 |    4.640e-13  |      1.1      |
    +--------+---------+---------+-----------+---------------+---------------+
    | 700.0  | 630.530 | 1507.61 | 1.1918e-9 |    1.537e-13  |      0.0      |
    +--------+---------+---------+-----------+---------------+---------------+

    """

    def __init__(self):

        super().__init__(TABLE_COESA62)

    def temperature(self, alt, geometric=True):
        """ Solves temperature at given altitude.

        Parameters
        ----------
        alt: ~astropy.units.Quantity
            Geometric/Geopotential altitude.
        geometric: boolean
            If `True`, assumes geometric altitude.

        Returns
        -------
        T: ~astropy.units.Quantity
            Kinetic temperature.
        """

        # Check if valid range and convert to geopotential
        Z, H = self._get_heights(alt, geometric=geometric)
        Zb, Hb, Tb, Lb, _ = self._get_base_parameters(H)

        if Z <= 90 * u.km:
            T = Tb + Lb * (H - Hb)
        else:
            T = Tb + Lb * (Z - Zb)

        return T.to(u.K)

    def pressure(self, alt, geometric=True):
        """ Solves pressure at given altitude.

        Parameters
        ----------
        alt: ~astropy.units.Quantity
            Geometric/Geopotential altitude.
        geometric: boolean
            If `True`, assumes geometric altitude.

        Returns
        -------
        p: ~astropy.units.Quantity
            Pressure at given altitude.
        """

        # Check if valid range and convert to geopotential
        Z, H = self._get_heights(alt, geometric=geometric)
        Zb, Hb, Tb, Lb, pb = self._get_base_parameters(H)
        g = gravity(Z, R_earth)

        # If Z > 90km, different formulas apply
        if Z <= 90 * u.km:
            if Lb == 0.0:
                p = pb * np.exp(-g * (H - Hb) / Tb / R_air)
            else:
                T = self.temperature(H, geometric=False)
                p = pb * (T / Tb) ** (-g / R_air / Lb)
        else:
            # TODO: Equation (1.2.10) should be applied avobe 90km
            raise NotImplementedError(
                "Pressure in COESA62 has just been implemented up to 90km."
            )

        return p

    def density(self, alt, geometric=True):
        """ Solves density at given altitude.

        Parameters
        ----------
        alt: ~astropy.units.Quantity
            Geometric/Geopotential altitude.
        geometric: boolean
            If `True`, assumes geometric altitude.

        Returns
        -------
        rho: ~astropy.units.Quantity
            Density at given altitude.
        """

        # Check if valid range
        Z, H = self._get_heights(alt, geometric=geometric)

        # TODO: implement atmosphere up to 1000km
        if Z > 90 * u.km:
            raise NotImplementedError(
                "Density in COESA62 has just been implemented up to 90km."
            )

        # Solve temperature and pressure
        T = self.temperature(H, geometric=False)
        p = self.pressure(H, geometric=False)
        rho = p / R_air / T
        return rho.to(u.kg / u.m ** 3)

    def properties(self, alt, geometric=True):
        """ Solves density at given height.

        Parameters
        ----------
        alt: ~astropy.units.Quantity
            Geometric/Geopotential height.
        geometric: boolean
            If `True`, assumes that `alt` argument is geometric kind.

        Returns
        -------
        T: ~astropy.units.Quantity
            Temperature at given height.
        p: ~astropy.units.Quantity
            Pressure at given height.
        rho: ~astropy.units.Quantity
            Density at given height.
        """
        T = self.temperature(alt, geometric=geometric)
        p = self.pressure(alt, geometric=geometric)
        rho = self.density(alt, geometric=geometric)

        return T, p, rho


# COESA76 TABLE
# Lower layers (Z < 86km)
Hb_lower76 = [0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0] * u.km
Zb_lower76 = [
    geopotential_to_geometric(H, R_earth).to(u.km).value for H in Hb_lower76
] * u.km

# Upper layers (Z >= 86km)
Zb_upper76 = [86.0, 91.0, 110.0, 120.0, 500.0, 1000.0] * u.km
Hb_upper76 = [
    geometric_to_geopotential(Z, R_earth).to(u.km).value for Z in Zb_upper76
] * u.km

# Geopotential and geometric layers
Hb_table76 = np.append(Hb_lower76, Hb_upper76).value * u.km
Zb_table76 = np.append(Zb_lower76, Zb_upper76).value * u.km

# Temperature layer
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

# Temperature gradient
Lb_table76 = (
    [-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0, 0.0, 12.0, 0.0, 0.0, 0.0] * u.K / u.km
)

# Pressure layer
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

# Generate table
TABLE_COESA76 = [Zb_table76, Hb_table76, Tb_table76, Lb_table76, pb_table76]

# Pressure and density coefficients
Z_coeff = [86, 91, 100, 110, 120, 150, 200, 300, 500, 750, 1000] * u.km
A_coeff_p = [
    0,
    0,
    0,
    0,
    2.283506e-07,
    1.209434e-08,
    8.113942e-10,
    9.814674e-11,
    -7.835161e-11,
    2.813255e-11,
    2.813255e-11,
]
B_coeff_p = [
    2.159582e-06,
    3.304895e-05,
    6.693926e-05,
    -6.539316e-05,
    -1.343221e-04,
    -9.692458e-06,
    -9.822568e-07,
    -1.654439e-07,
    1.964589e-07,
    -1.120689e-07,
    -1.120689e-07,
]
C_coeff_p = [
    -4.836957e-04,
    -0.009062730,
    -0.01945388,
    0.02485568,
    0.02999016,
    0.003002041,
    4.687616e-04,
    1.148115e-04,
    -1.657213e-04,
    1.695568e-04,
    1.695568e-04,
]
D_coeff_p = [
    -0.1425192,
    0.6516698,
    1.719080,
    -3.223620,
    -3.055446,
    -0.4523015,
    -0.1231710,
    -0.05431334,
    0.04305869,
    -0.1188941,
    -0.1188941,
]
E_coeff_p = [
    13.47530,
    -11.03037,
    -47.75030,
    135.9355,
    113.5764,
    19.19151,
    3.067409,
    -2.011365,
    -14.77132,
    14.56718,
    14.56718,
]

A_coeff_rho = [
    0,
    0,
    -1.240774e-05,
    0,
    3.661771e-07,
    1.906032e-08,
    1.199282e-09,
    1.140564e-10,
    8.105631e-12,
    -3.701195e-12,
    -3.701195e-12,
]

B_coeff_rho = [
    -3.322622e-06,
    2.873405e-05,
    0.005162063,
    -8.854164e-05,
    -2.154344e-04,
    -1.527799e-05,
    -1.451051e-06,
    -2.130756e-07,
    -2.358417e-09,
    -8.608611e-09,
    -8.608611e-09,
]

C_coeff_rho = [
    9.111460e-04,
    -0.008492037,
    -0.8048342,
    0.03373254,
    0.04809214,
    0.004724294,
    6.910474e-04,
    1.570762e-04,
    -2.635110e-06,
    5.118829e-05,
    5.118829e-05,
]

D_coeff_rho = [
    -0.2609971,
    0.6541179,
    55.55996,
    -4.390837,
    -4.884744,
    -0.6992340,
    -0.1736220,
    -0.07029296,
    -0.01562608,
    -0.06600998,
    -0.06600998,
]

E_coeff_rho = [
    5.944694,
    -23.62010,
    -1443.338,
    176.5294,
    172.3597,
    20.50921,
    -5.321644,
    -12.89844,
    -20.02246,
    -6.137674,
    -6.137674,
]

p_coeff = [A_coeff_p, B_coeff_p, C_coeff_p, D_coeff_p, E_coeff_p]
rho_coeff = [A_coeff_rho, B_coeff_rho, C_coeff_rho, D_coeff_rho, E_coeff_rho]


class COESA76(COESA):
    """ U.S Standard Atmosphere 1976.

    +--------+---------+---------+-----------+---------------+---------------+
    | Z (km) |  H (km) |  T (K)  |  p (mbar) | rho (kg / m3) | beta (K / km) |
    +--------+---------+---------+-----------+---------------+---------------+
    |   0.0  |   0.0   | 288.150 | 1.01325e3 |     1.2250    |      -6.5     |
    +--------+---------+---------+-----------+---------------+---------------+
    | 11.019 |   11.0  | 216.650 |  2.2632e2 |   3.6392e-1   |      0.0      |
    +--------+---------+---------+-----------+---------------+---------------+
    | 20.063 |   20.0  | 216.650 |  5.4748e1 |   8.8035e-2   |      1.0      |
    +--------+---------+---------+-----------+---------------+---------------+
    | 32.162 |   32.0  | 288.650 |  8.6801e0 |   1.3225e-2   |      2.8      |
    +--------+---------+---------+-----------+---------------+---------------+
    | 47.350 |   47.0  | 270.650 |  1.1090e0 |   1.4275e-3   |      0.0      |
    +--------+---------+---------+-----------+---------------+---------------+
    | 51.413 |   51.0  | 270.650 | 6.6938e-1 |   8.6160e-4   |      -2.8     |
    +--------+---------+---------+-----------+---------------+---------------+
    | 71.802 |   71.0  | 214.650 | 3.9564e-2 |   6.4211e-5   |      -2.0     |
    +--------+---------+---------+-----------+---------------+---------------+
    |  86.0  | 84.8520 |  186.87 | 3.7338e-3 |    6.958e-6   |      0.0      |
    +--------+---------+---------+-----------+---------------+---------------+
    |  91.0  |  89.716 |  186.87 | 1.5381e-3 |    2.860e-6   |   elliptical  |
    +--------+---------+---------+-----------+---------------+---------------+
    |  110.0 | 108.129 |  240.00 | 7.1042e-5 |    9.708e-8   |      12.0     |
    +--------+---------+---------+-----------+---------------+---------------+
    |  120.0 | 117.777 |  360.00 | 2.5382e-5 |    2.222e-8   |  exponential  |
    +--------+---------+---------+-----------+---------------+---------------+
    |  500.0 | 463.540 |  999.24 | 3.0236e-9 |   5.215e-13   |  exponential  |
    +--------+---------+---------+-----------+---------------+---------------+
    | 1000.0 | 864.071 |   1000  | 7.5138e-5 |   3.561e-15   |  exponential  |
    +--------+---------+---------+-----------+---------------+---------------+

    """

    def __init__(self):

        super().__init__(TABLE_COESA76)

    def _get_coefficients_avobe_86(self, Z, table_coeff):
        """ Returns corresponding coefficients for 4th order polynomial approximation.

        Parameters
        ----------
        Z: ~astropy.units.Quantity
            Geometric altitude
        table_coeff: list
            List containing different coefficient lists.

        Returns
        -------
        coeff_list: list
            List of corresponding coefficients
        """
        # Get corresponding coefficients
        i = self._get_index(Z, Z_coeff)
        coeff_list = []
        for coeff_set in table_coeff:
            x = coeff_set[i]
            coeff_list.append(x)

        return coeff_list

    def temperature(self, alt, geometric=True):
        """ Returns temperature at given height.

        Parameters
        ----------
        alt: ~astropy.units.Quantity
            Geometric/Geopotential height.
        geometric: boolean
            If `True`, assumes that `alt` argument is geometric kind.

        Returns
        -------
        T: ~astropy.units.Quantity
            Temperature at given height.
        """
        # Check if valid range and convert to geopotential
        Z, H = self._get_heights(alt, geometric=geometric)
        Zb, Hb, Tb, Lb, _ = self._get_base_parameters(H)

        # Different altitudes imply different equations
        # TODO: Include mean air molecular weight corrections between [80, 86]km
        if Z <= self._Zb_table[7]:
            # Equation (23)
            T = Tb + Lb * (H - Hb)
        elif Z > self._Zb_table[7] and Z <= self._Zb_table[8]:
            # Equation (25)
            T = 186.8673 * u.K
        elif Z > self._Zb_table[8] and Z <= self._Zb_table[9]:
            # Equation (27)
            Tc = 263.1905 * u.K
            A = -76.3232 * u.K
            a = -19.9429 * u.km
            T = Tc + A * (1 - ((Z - self._Zb_table[8]) / a) ** 2) ** 0.5
        elif Z >= self._Zb_table[9] and Z <= self._Zb_table[10]:
            # Equation (29)
            T = Tb + Lb * (Z - Zb)
        else:
            # Equation (32)
            _lambda = 0.01875 / u.km
            epsilon = (
                (Z - self._Zb_table[10])
                * (R_earth + self._Zb_table[10])
                / (R_earth + Z)
            )
            Tinf = 1000 * u.K
            T = Tinf - (Tinf - self._Tb_table[10]) * np.exp(-_lambda * epsilon)

        return T.to(u.K)

    def pressure(self, alt, geometric=True):
        """ Returns pressure at given height.

        Parameters
        ----------
        alt: ~astropy.units.Quantity
            Geometric/Geopotential height.
        geometric: boolean
            If `True`, assumes that `alt` argument is geometric kind.

        Returns
        -------
        p: ~astropy.units.Quantity
            Pressure at given height.
        """
        # Checks are done in temperature function
        Z, H = self._get_heights(alt, geometric=geometric)

        T = self.temperature(H, geometric=False)
        Zb, Hb, Tb, Lb, pb = self._get_base_parameters(H)

        # Obtain gravity magnitude
        g = gravity(Z, R_earth)

        # If above 86[km] usual formulation is applied
        if Z <= 86 * u.km:
            if Lb == 0.0:
                p = pb * np.exp(-g * (H - Hb) / Tb / R_air)
            else:
                p = pb * (T / Tb) ** (-g / R_air / Lb)
            return p
        else:
            # A 4th order polynomial is used to approximate pressure.  This was
            # directly taken from: http://www.braeunig.us/space/atmmodel.htm
            A, B, C, D, E = self._get_coefficients_avobe_86(Z, p_coeff)

            # Solve the polynomial
            z = Z.to(u.km).value
            return np.exp(A * z ** 4 + B * z ** 3 + C * z ** 2 + D * z + E) * u.Pa

    def density(self, alt, geometric=True):
        """ Solves density at given height.

        Parameters
        ----------
        alt: ~astropy.units.Quantity
            Geometric/Geopotential height.
        geometric: boolean
            If `True`, assumes that `alt` argument is geometric kind.

        Returns
        -------
        rho: ~astropy.units.Quantity
            Density at given height.
        """
        # Check if valid range
        Z, H = self._get_heights(alt, geometric=geometric)

        # Solve temperature and pressure
        if alt <= 86 * u.km:
            T = self.temperature(H, geometric=False)
            p = self.pressure(H, geometric=False)
            rho = p / R_air / T
            return rho.to(u.kg / u.m ** 3)
        else:
            # A 4th order polynomial is used to approximate pressure.  This was
            # directly taken from: http://www.braeunig.us/space/atmmodel.htm
            A, B, C, D, E = self._get_coefficients_avobe_86(Z, rho_coeff)

            # Solve the polynomial
            z = Z.to(u.km).value
            return (
                np.exp(A * z ** 4 + B * z ** 3 + C * z ** 2 + D * z + E)
                * u.kg
                / u.m ** 3
            )

    def properties(self, alt, geometric=True):
        """ Solves density at given height.

        Parameters
        ----------
        alt: ~astropy.units.Quantity
            Geometric/Geopotential height.
        geometric: boolean
            If `True`, assumes that `alt` argument is geometric kind.

        Returns
        -------
        T: ~astropy.units.Quantity
            Temperature at given height.
        p: ~astropy.units.Quantity
            Pressure at given height.
        rho: ~astropy.units.Quantity
            Density at given height.
        """
        T = self.temperature(alt, geometric=geometric)
        p = self.pressure(alt, geometric=geometric)
        rho = self.density(alt, geometric=geometric)
        return T, p, rho
