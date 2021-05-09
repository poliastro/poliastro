""" The U.S. Standard Atmosphere 1976 is an idealized, steady-state model of
mean annual conditions of Earth's atmosphere from the surface to 1000 km at
latitude 45N, as it is assumed to exist during a period with moderate solar
activity. The defining meteorological elements are sea-level temperature and
pressure, and a temperature-height profile to 1000 km. The air is assumed to be
dry, and at heights sufficiently below 86 km, the atmosphere is assumed to be
homogeneously mixed with a relative-volume composition leading to a constant
mean molecular weight.

Since 1976 many constants such us Earth's radius or Avogadro's number have been
updated. In order to have a pure COESA76 atmospheric model, the official paper
values were used.

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

import numpy as np
from astropy import units as u
from astropy.io import ascii
from astropy.utils.data import get_pkg_data_filename

from poliastro.earth.atmosphere.base import COESA

# Following constants come from original U.S Atmosphere 1962 paper so a pure
# model of this atmosphere can be implemented
R = 8314.32 * u.J / u.kmol / u.K
R_air = 287.053 * u.J / u.kg / u.K
k = 1.380622e-23 * u.J / u.K
Na = 6.022169e-26 / u.kmol
g0 = 9.80665 * u.m / u.s ** 2
r0 = 6356.766 * u.km
M0 = 28.9644 * u.kg / u.kmol
P0 = 101325 * u.Pa
T0 = 288.15 * u.K
Tinf = 1000 * u.K
gamma = 1.4
alpha = 34.1632 * u.K / u.km
beta = 1.458e-6 * (u.kg / u.s / u.m / (u.K) ** 0.5)
S = 110.4 * u.K

# Reading layer parameters file
coesa76_data = ascii.read(get_pkg_data_filename("data/coesa76.dat"))
b_levels = coesa76_data["b"].data
zb_levels = coesa76_data["Zb [km]"].data * u.km
hb_levels = coesa76_data["Hb [km]"].data * u.km
Tb_levels = coesa76_data["Tb [K]"].data * u.K
Lb_levels = coesa76_data["Lb [K/km]"].data * u.K / u.km
pb_levels = coesa76_data["pb [mbar]"].data * u.mbar

# Reading pressure and density coefficients files
p_data = ascii.read(get_pkg_data_filename("data/coesa76_p.dat"))
rho_data = ascii.read(get_pkg_data_filename("data/coesa76_rho.dat"))

# Zip coefficients for each altitude
z_coeff = p_data["z [km]"].data * u.km
p_coeff = [
    p_data["A"].data,
    p_data["B"].data,
    p_data["C"].data,
    p_data["D"].data,
    p_data["E"].data,
]
rho_coeff = [
    rho_data["A"].data,
    rho_data["B"].data,
    rho_data["C"].data,
    rho_data["D"].data,
    rho_data["E"].data,
]


class COESA76(COESA):
    """Holds the model for U.S Standard Atmosphere 1976."""

    def __init__(self):
        """Constructor for the class."""
        super().__init__(
            b_levels, zb_levels, hb_levels, Tb_levels, Lb_levels, pb_levels
        )

    def _get_coefficients_avobe_86(self, z, table_coeff):
        """Returns corresponding coefficients for 4th order polynomial approximation.

        Parameters
        ----------
        z: ~astropy.units.Quantity
            Geometric altitude
        table_coeff: list
            List containing different coefficient lists.

        Returns
        -------
        coeff_list: list
            List of corresponding coefficients
        """

        # Get corresponding coefficients
        i = self._get_index(z, z_coeff)
        coeff_list = []
        for X_set in table_coeff:
            coeff_list.append(X_set[i])

        return coeff_list

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
        Tb = self.Tb_levels[i]
        Lb = self.Lb_levels[i]
        hb = self.hb_levels[i]

        # Apply different equations
        if z < self.zb_levels[7]:
            # Below 86km
            # TODO: Apply air mean molecular weight ratio factor
            Tm = Tb + Lb * (h - hb)
            T = Tm
        elif z >= self.zb_levels[7] and z < self.zb_levels[8]:
            # [86km, 91km)
            T = 186.87 * u.K
        elif z >= self.zb_levels[8] and z < self.zb_levels[9]:
            # [91km, 110km]
            Tc = 263.1905 * u.K
            A = -76.3232 * u.K
            a = -19.9429 * u.km
            T = Tc + A * (1 - ((z - self.zb_levels[8]) / a) ** 2) ** 0.5
        elif z >= self.zb_levels[9] and z < self.zb_levels[10]:
            # [110km, 120km]
            T = 240 * u.K + Lb * (z - self.zb_levels[9])
        else:
            T10 = 360.0 * u.K
            _gamma = self.Lb_levels[9] / (Tinf - T10)
            epsilon = (z - self.zb_levels[10]) * (r0 + self.zb_levels[10]) / (r0 + z)
            T = Tinf - (Tinf - T10) * np.exp(-_gamma * epsilon)

        return T.to(u.K)

    def pressure(self, alt, geometric=True):
        """Solves pressure at given altitude.

        Parameters
        ----------
        alt: ~astropy.units.Quantity
            Geometric/Geopotential altitude.
        geometric: bool
            If `True`, assumes geometric altitude kind.

        Returns
        -------
        p: ~astropy.units.Quantity
            Pressure at given altitude.
        """
        # Test if altitude is inside valid range
        z, h = self._check_altitude(alt, r0, geometric=geometric)

        # Obtain gravity magnitude
        # Get base parameters
        i = self._get_index(z, self.zb_levels)
        Tb = self.Tb_levels[i]
        Lb = self.Lb_levels[i]
        hb = self.hb_levels[i]
        pb = self.pb_levels[i]

        # If above 86[km] usual formulation is applied
        if z < 86 * u.km:
            if Lb == 0.0 * u.K / u.km:
                p = pb * np.exp(-alpha * (h - hb) / Tb)
            else:
                T = self.temperature(z)
                p = pb * (Tb / T) ** (alpha / Lb)
        else:
            # TODO: equation (33c) should be applied instead of using coefficients

            # A 4th order polynomial is used to approximate pressure.  This was
            # directly taken from: http://www.braeunig.us/space/atmmodel.htm
            A, B, C, D, E = self._get_coefficients_avobe_86(z, p_coeff)

            # Solve the polynomial
            z = z.to(u.km).value
            p = np.exp(A * z ** 4 + B * z ** 3 + C * z ** 2 + D * z + E) * u.Pa

        return p.to(u.Pa)

    def density(self, alt, geometric=True):
        """Solves density at given height.

        Parameters
        ----------
        alt: ~astropy.units.Quantity
            Geometric/Geopotential height.
        geometric: bool
            If `True`, assumes that `alt` argument is geometric kind.

        Returns
        -------
        rho: ~astropy.units.Quantity
            Density at given height.
        """
        # Test if altitude is inside valid range
        z, h = self._check_altitude(alt, r0, geometric=geometric)

        # Solve temperature and pressure
        if z <= 86 * u.km:
            T = self.temperature(z)
            p = self.pressure(z)
            rho = p / R_air / T
        else:
            # TODO: equation (42) should be applied instead of using coefficients

            # A 4th order polynomial is used to approximate pressure.  This was
            # directly taken from: http://www.braeunig.us/space/atmmodel.htm
            A, B, C, D, E = self._get_coefficients_avobe_86(z, rho_coeff)

            # Solve the polynomial
            z = z.to(u.km).value
            rho = (
                np.exp(A * z ** 4 + B * z ** 3 + C * z ** 2 + D * z + E)
                * u.kg
                / u.m ** 3
            )

        return rho.to(u.kg / u.m ** 3)

    def properties(self, alt, geometric=True):
        """Solves temperature, pressure, density at given height.

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

        if z > 86 * u.km:
            raise ValueError(
                "Speed of sound in COESA76 has just been implemented up to 86km."
            )
        T = self.temperature(alt, geometric).value
        # Using eqn-(50)
        Cs = ((gamma * R_air.value * T) ** 0.5) * (u.m / u.s)

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

        if z > 86 * u.km:
            raise ValueError(
                "Dynamic Viscosity in COESA76 has just been implemented up to 86km."
            )
        T = self.temperature(alt, geometric).value
        # Using eqn-(51)
        mu = (beta.value * T ** 1.5 / (T + S.value)) * (u.N * u.s / (u.m) ** 2)

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

        if z > 86 * u.km:
            raise ValueError(
                "Thermal conductivity in COESA76 has just been implemented up to 86km."
            )
        T = self.temperature(alt, geometric=geometric).value
        # Using eqn-(53)
        k = (2.64638e-3 * T ** 1.5 / (T + 245.4 * (10 ** (-12.0 / T)))) * (
            u.J / u.m / u.s / u.K
        )

        return k
