""" The U.S. Standard Atmosphere 1966 depicts idealized middle-latitude
year-round mean conditions for the range of solar activity that occurs between
sunspot minimum and sunspot maximum.

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

import numpy as np
from astropy import units as u
from astropy.io import ascii
from astropy.units import imperial
from astropy.utils.data import get_pkg_data_filename
from scipy.integrate import quad

from poliastro.earth.atmosphere.base import COESA

# Constants come from the original paper to achieve pure implementation
r0 = 6356.766 * u.km
p0 = 1.013250e5 * u.Pa
rho0 = 1.2250 * u.K
T0 = 288.15 * u.K
g0 = 9.80665 * u.m / u.s ** 2
S = 110.4 * u.K
Ti = 273.15 * u.K
beta = 1.458e-6 * u.kg / u.s / u.m / u.K ** (0.5)
_gamma = 1.4
sigma = 3.65e-10 * u.m
N = 6.02257e26 * (u.kg * u.mol) ** -1
R = 8314.32 * u.J / u.kmol / u.K
R_air = 287.053 * u.J / u.kg / u.K
alpha = 34.1632 * u.K / u.km

# Reading layer parameters file
coesa_file = get_pkg_data_filename("data/coesa62.dat")
coesa62_data = ascii.read(coesa_file)
b_levels = coesa62_data["b"].data
zb_levels = coesa62_data["Zb [km]"].data * u.km
hb_levels = coesa62_data["Hb [km]"].data * u.km
Tb_levels = coesa62_data["Tb [K]"].data * u.K
Lb_levels = coesa62_data["Lb [K/km]"].data * u.K / u.km
pb_levels = coesa62_data["pb [mbar]"].data * u.mbar


class COESA62(COESA):
    """Holds the model for U.S Standard Atmosphere 1962."""

    def __init__(self):
        """Constructor for the class."""
        super().__init__(
            b_levels, zb_levels, hb_levels, Tb_levels, Lb_levels, pb_levels
        )

    def temperature(self, alt, geometric=True):
        """Solves for temperature at given altitude.

        Parameters
        ----------
        alt: ~astropy.units.Quantity
            Geometric/Geopotential altitude.
        geometric: bool
            If `True`, assumes geometric altitude kind.

        Returns
        -------
        T: ~astropy.units.Quantity
            Kinetic temeperature.
        """

        # Test if altitude is inside valid range
        z, h = self._check_altitude(alt, r0, geometric=geometric)

        # Get base parameters
        i = self._get_index(z, self.zb_levels)
        zb = self.zb_levels[i]
        Tb = self.Tb_levels[i]
        Lb = self.Lb_levels[i]
        hb = self.hb_levels[i]

        # Apply different equations
        if z <= 90 * u.km:
            T = Tb + Lb * (h - hb)
        else:
            T = Tb + Lb * (z - zb)

        return T.to(u.K)

    def pressure(self, alt, geometric=True):
        """Solves pressure at given altitude.

        Parameters
        ----------
        alt: ~astropy.units.Quantity
            Geometric/Geopotential altitude.
        geometric: bool
            If `True`, assumes geometric altitude.

        Returns
        -------
        p: ~astropy.units.Quantity
            Pressure at given altitude.
        """

        # Check if valid range and convert to geopotential
        z, h = self._check_altitude(alt, r0, geometric=geometric)

        # Get base parameters
        i = self._get_index(z, self.zb_levels)
        zb = self.zb_levels[i]
        hb = self.hb_levels[i]
        Tb = self.Tb_levels[i]
        Lb = self.Lb_levels[i]
        pb = self.pb_levels[i]

        # If z <= 90km then apply eqn 1.2.10-(3)
        if z <= 90 * u.km:
            # If Lb is zero then apply eqn 1.2.10-(4)
            if Lb == 0.0:
                p = pb * np.exp(-g0 * (h - hb) / Tb / R_air)
            else:
                T = self.temperature(z)
                p = pb * (T / Tb) ** (-g0 / R_air / Lb)

        # If 90 < Z < 700 km then eqn 1.2.10-(5) is applied
        else:
            # Converting all the units into SI unit and taking their magnitude
            Lb_v = Lb.to(u.K / u.m).value
            r0_v = r0.to(u.m).value
            z_v = z.to(u.m).value
            zb_v = zb.to(u.m).value
            Tb_v = Tb.value
            g0_v = g0.value
            R_air_v = R_air.value

            # Putting g = (g0*(r0/(r0 +z))**2) in (g * dz / z - zb + Tb/Lb)
            # and integrating it.
            integrand = quad(
                lambda x: (g0_v * (r0_v / (r0_v + x)) ** 2) / (x - zb_v + Tb_v / Lb_v),
                zb_v,
                z_v,
            )

            pb = pb.to(u.Pa)
            p = (pb * np.exp((-1 / R_air_v / Lb_v) * integrand[0])).to(u.mbar)

        return p

    def density(self, alt, geometric=True):
        """Solves density at given altitude.

        Parameters
        ----------
        alt: ~astropy.units.Quantity
            Geometric/Geopotential altitude.
        geometric: bool
            If `True`, assumes geometric altitude.

        Returns
        -------
        rho: ~astropy.units.Quantity
            Density at given altitude.
        """

        # Check if valid range and convert to geopotential
        z, h = self._check_altitude(alt, r0, geometric=geometric)

        # Solve temperature and pressure
        T = self.temperature(z)
        p = self.pressure(z)
        rho = p / R_air / T

        return rho.to(u.kg / u.m ** 3)

    def properties(self, alt, geometric=True):
        """Solves density at given height.

        Parameters
        ----------
        alt: ~astropy.units.Quantity
            Geometric/Geopotential height.
        geometric: bool
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

    def sound_speed(self, alt, geometric=True):
        """Solves speed of sound at given height.

        Parameters
        ----------
        alt: ~astropy.units.Quantity
            Geometric/Geopotential height.
        geometric: bool
            If `True`, assumes that `alt` argument is geometric kind.

        Returns
        -------
        Cs: ~astropy.units.Quantity
            Speed of Sound at given height.
        """
        # Check if valid range and convert to geopotential
        z, h = self._check_altitude(alt, r0, geometric=geometric)

        if z > 90 * u.km:
            raise ValueError(
                "Speed of sound in COESA62 has just been implemented up to 90km."
            )
        T = self.temperature(alt, geometric).value
        # Using eqn-1.3.7-(1)
        Cs = ((_gamma * R_air.value * T) ** 0.5) * (u.m / u.s)

        return Cs

    def viscosity(self, alt, geometric=True):
        """Solves dynamic viscosity at given height.

        Parameters
        ----------
        alt: ~astropy.units.Quantity
            Geometric/Geopotential height.
        geometric: bool
            If `True`, assumes that `alt` argument is geometric kind.

        Returns
        -------
        mu: ~astropy.units.Quantity
            Dynamic viscosity at given height.
        """
        # Check if valid range and convert to geopotential
        z, h = self._check_altitude(alt, r0, geometric=geometric)

        if z > 90 * u.km:
            raise ValueError(
                "Dynamic Viscosity in COESA62 has just been implemented up to 90km."
            )
        T = self.temperature(alt, geometric).value
        # Using eqn-1.3.8-(1)
        mu = (beta.value * T ** 1.5 / (T + S.value)) * (u.kg / u.m / u.s)

        return mu

    def thermal_conductivity(self, alt, geometric=True):
        """Solves coefficient of thermal conductivity at given height.

        Parameters
        ----------
        alt: ~astropy.units.Quantity
            Geometric/Geopotential height.
        geometric: bool
            If `True`, assumes that `alt` argument is geometric kind.

        Returns
        -------
        k: ~astropy.units.Quantity
            coefficient of thermal conductivity at given height.
        """
        # Check if valid range and convert to geopotential
        z, h = self._check_altitude(alt, r0, geometric=geometric)

        if z > 90 * u.km:
            raise ValueError(
                "Thermal conductivity in COESA62 has just been implemented up to 90km."
            )

        T = self.temperature(alt, geometric=geometric).value
        # Using eqn-1.3.10-(1)
        k = (6.325e-7 * T ** 1.5 / (T + 245.4 * (10 ** (-12.0 / T)))) * (
            imperial.kcal / u.m / u.s / u.K
        )

        return k
