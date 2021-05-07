"""
Given an exospheric temperature, Jacchia77 returns model
atmospheric altitude profiles of temperature, the number
densities of N2, O2, O, Ar, He, H, the sum thereof, and the
molecular weight.

For altitudes of 90 km and above, we use the 1977 model of
Jacchia [Ja77].  H-atom densities are returned as non-zero
for altitudes of 150 km and above if the maximum altitude
requested is 500 km or more.

REFERENCES:

Ja77    L. G. Jacchia, "Thermospheric Temperature, Density
        and Composition: New Models," SAO Special Report No.
        375 (Smithsonian Institution Astrophysical
        Observatory, Cambridge, MA, March 15, 1977).

Fortran Implementation:
https://ccmc.gsfc.nasa.gov/pub/modelweb/atmospheric/jacchia/jacchia-77/

"""
import numpy as np
from astropy import units as u

# Following constants have been taken from the fortran implementation
pi2 = np.pi / 2
wm0 = 28.96
wmN2 = 28.0134
wmO2 = 31.9988
wmO = 15.9994
wmAr = 39.948
wmHe = 4.0026
wmH = 1.0079
qN2 = 0.78110
qO2 = 0.20955
qAr = 0.009343
qHe = 0.000005242
R0 = 6356.766
R = 8314.32 * u.J / (u.kg * u.mol)
k = 1.380622e-23 * u.J / u.K


class Jacchia77:
    """Holds the model for U.S Standard Atmosphere 1962."""

    def __init__(self):
        self.E5M = [0.0 for _ in range(11)]
        self.E6P = [0.0 for _ in range(11)]
        self.x = 0.0
        self.y = 0.0

    def _altitude_profile(self, alt, Texo):
        # Raise Value Error if alt < 90 km or alt > 2500 km.
        if alt < 90 or alt > 2500:
            raise ValueError("Jacchia77 has been implemented in range 90km - 2500km.")

        alt = int(
            alt + 1
        )  # in fortran the upper limits are included. in python are not.
        Texo = int(Texo)

        self.Z = [0.0 for _ in range(alt)]
        self.T = [0.0 for _ in range(alt)]
        self.CN2 = [0.0 for _ in range(alt)]
        self.CO2 = [0.0 for _ in range(alt)]
        self.CO = [0.0 for _ in range(alt)]
        self.CAr = [0.0 for _ in range(alt)]
        self.CHe = [0.0 for _ in range(alt)]
        self.CH = [0.0 for _ in range(alt)]
        self.CM = [0.0 for _ in range(alt)]
        self.WM = [0.0 for _ in range(alt)]

        for iz in range(90, alt):
            self.Z[iz] = iz
            self.CH[iz] = 0

            if iz <= 90:
                self.T[iz] = 188
            elif Texo < 188.1:
                self.T[iz] = 188
            else:
                self.x = 0.0045 * (Texo - 188.0)
                Tx = 188 + 110.5 * np.log(self.x + np.sqrt(self.x * self.x + 1))
                Gx = pi2 * 1.9 * (Tx - 188.0) / (125.0 - 90.0)
                if iz <= 125:
                    self.T[iz] = Tx + ((Tx - 188.0) / pi2) * np.arctan(
                        (Gx / (Tx - 188.0))
                        * (self.Z[iz] - 125.0)
                        * (
                            1.0
                            + 1.7 * ((self.Z[iz] - 125.0) / (self.Z[iz] - 90.0)) ** 2
                        )
                    )
                else:
                    self.T[iz] = Tx + ((Texo - Tx) / pi2) * np.arctan(
                        (Gx / (Texo - Tx))
                        * (self.Z[iz] - 125.0)
                        * (1.0 + 5.5e-5 * (self.Z[iz] - 125.0) ** 2)
                    )
            if iz <= 100:
                self.x = iz - 90
                self.E5M[iz - 90] = 28.89122 + self.x * (
                    -2.83071e-2
                    + self.x
                    * (
                        -6.59924e-3
                        + self.x
                        * (
                            -3.39574e-4
                            + self.x * (+6.19256e-5 + self.x * (-1.84796e-6))
                        )
                    )
                )
                if iz <= 90:
                    self.E6P[0] = 7.145e13 * self.T[90]
                else:
                    G0 = (1 + self.Z[iz - 1] / R0) ** (-2)
                    G1 = (1 + self.Z[iz] / R0) ** (-2)
                    self.E6P[iz - 90] = self.E6P[iz - 91] * np.exp(
                        -0.5897446
                        * (
                            G1 * self.E5M[iz - 90] / self.T[iz]
                            + G0 * self.E5M[iz - 91] / self.T[iz - 1]
                        )
                    )

                self.x = self.E5M[iz - 90] / wm0
                self.y = self.E6P[iz - 90] / self.T[iz]

                self.CN2[iz] = qN2 * self.y * self.x
                self.CO[iz] = 2.0 * (1.0 - self.x) * self.y
                self.CO2[iz] = (self.x * (1.0 + qO2) - 1.0) * self.y
                self.CAr[iz] = qAr * self.y * self.x
                self.CHe[iz] = qHe * self.y * self.x
                self.CH[iz] = 0
            else:
                G0 = (1 + self.Z[iz - 1] / R0) ** (-2)
                G1 = (1 + self.Z[iz] / R0) ** (-2)

                self.x = 0.5897446 * (G1 / self.T[iz] + G0 / self.T[iz - 1])
                self.y = self.T[iz - 1] / self.T[iz]
                self.CN2[iz] = self.CN2[iz - 1] * self.y * np.exp(-wmN2 * self.x)
                self.CO2[iz] = self.CO2[iz - 1] * self.y * np.exp(-wmO2 * self.x)
                self.CO[iz] = self.CO[iz - 1] * self.y * np.exp(-wmO * self.x)
                self.CAr[iz] = self.CAr[iz - 1] * self.y * np.exp(-wmAr * self.x)
                self.CHe[iz] = (
                    self.CHe[iz - 1] * (self.y ** 0.62) * np.exp(-wmHe * self.x)
                )
                self.CH[iz] = 0

        self._O_and_O2_correction(alt, Texo)

        if alt >= 500:
            self._H_correction(alt, Texo)

        return (
            self.Z * u.km,
            self.T * u.K,
            np.array(self.CN2) * 1e6 * (u.m) ** -3,
            np.array(self.CO2) * 1e6 * (u.m) ** -3,
            np.array(self.CO) * 1e6 * (u.m) ** -3,
            np.array(self.CAr) * 1e6 * (u.m) ** -3,
            np.array(self.CHe) * 1e6 * (u.m) ** -3,
            np.array(self.CH) * 1e6 * (u.m) ** -3,
            np.array(self.CM) * 1e6 * (u.m) ** -3,
            self.WM * (u.kg),
        )

    def _H_correction(self, alt, Texo):
        """Calculate [H] from Jacchia 1977 formulas"""

        phid00 = 10.0 ** (6.9 + 28.9 * Texo ** (-0.25)) / 2.0e20
        phid00 = phid00 * 5.24e2
        H_500 = 10.0 ** (-0.06 + 28.9 * Texo ** (-0.25))
        # print(alt)
        for iz in range(150, alt):
            phid0 = phid00 / np.sqrt(self.T[iz])
            self.WM[iz] = (
                wmH * 0.5897446 * ((1.0 + self.Z[iz] / R0) ** (-2)) / self.T[iz] + phid0
            )
            self.CM[iz] = self.CM[iz] * phid0

        self.y = self.WM[150]
        self.WM[150] = 0

        for iz in range(151, alt):
            self.x = self.WM[iz - 1] + (self.y + self.WM[iz])
            self.y = self.WM[iz]
            self.WM[iz] = self.x

        for iz in range(150, alt):
            self.WM[iz] = np.exp(self.WM[iz]) * (self.T[iz] / self.T[150]) ** 0.75
            self.CM[iz] = self.WM[iz] * self.CM[iz]

        self.y = self.CM[150]
        self.CM[150] = 0

        for iz in range(151, alt):
            self.x = self.CM[iz - 1] + 0.5 * (self.y + self.CM[iz])
            self.y = self.CM[iz]
            self.CM[iz] = self.x

        for iz in range(150, alt):
            self.CH[iz] = (self.WM[500] / self.WM[iz]) * (
                H_500 - (self.CM[iz] - self.CM[500])
            )

        for iz in range(150, alt):
            self.CM[iz] = (
                self.CN2[iz]
                + self.CO2[iz]
                + self.CO[iz]
                + self.CAr[iz]
                + self.CHe[iz]
                + self.CH[iz]
            )
            self.WM[iz] = (
                wmN2 * self.CN2[iz]
                + wmO2 * self.CO2[iz]
                + wmO * self.CO[iz]
                + wmAr * self.CAr[iz]
                + wmHe * self.CHe[iz]
                + wmH * self.CH[iz]
            ) / self.CM[iz]

    def _O_and_O2_correction(self, alt, Texo):
        """Add Jacchia 1977 empirical corrections to [O] and [O2]"""

        for iz in range(90, alt):
            self.CO2[iz] = self.CO2[iz] * (
                10.0 ** (-0.07 * (1.0 + np.tanh(0.18 * (self.Z[iz] - 111.0))))
            )
            self.CO[iz] = self.CO[iz] * (
                10.0 ** (-0.24 * np.exp(-0.009 * (self.Z[iz] - 97.7) ** 2))
            )
            self.CM[iz] = (
                self.CN2[iz]
                + self.CO2[iz]
                + self.CO[iz]
                + self.CAr[iz]
                + self.CHe[iz]
                + self.CH[iz]
            )
            self.WM[iz] = (
                wmN2 * self.CN2[iz]
                + wmO2 * self.CO2[iz]
                + wmO * self.CO[iz]
                + wmAr * self.CAr[iz]
                + wmHe * self.CHe[iz]
                + wmH * self.CH[iz]
            ) / self.CM[iz]

    def altitude_profile(self, alt, Texo):
        """Solves for atmospheric altitude profile at given altitude and exospheric temperature.

        Parameters
        ----------
        alt: ~astropy.units.Quantity
            Geometric/Geopotential altitude.
        Texo: ~astropy.units.Quantity
            Exospheric temperature

        Returns
        -------
        altitude_profile: list
            [altitude(Z), T, N2, O2, O, Ar, He, H, Total number density, Mean Molecular weight]
        """
        # checking if the units entered are km
        if alt.unit == u.km:
            if 150 <= alt.value < 500:
                alt_properties = self._altitude_profile(500, Texo.value)
            else:
                alt_properties = self._altitude_profile(alt.value, Texo.value)

        return [last[int(alt.value)] for last in alt_properties]

    def temperature(self, alt, Texo):
        """Solves for temperature at given altitude and exospheric temperature.

        Parameters
        ----------
        alt: ~astropy.units.Quantity
            Geometric/Geopotential altitude.
        Texo: ~astropy.units.Quantity
            Exospheric temperature

        Returns
        -------
        T: ~astropy.units.Quantity
            Absolute temeperature  and exospheric temperature
        """
        T = self.altitude_profile(alt, Texo)[1]
        return T

    def pressure(self, alt, Texo):
        """Solves pressure at given altitude and exospheric temperature.

        Parameters
        ----------
        alt: ~astropy.units.Quantity
            Geometric/Geopotential altitude.
        Texo: ~astropy.units.Quantity
            Exospheric temperature

        Returns
        -------
        p: ~astropy.units.Quantity
            Pressure at given altitude  and exospheric temperature.
        """
        alt_profile = self.altitude_profile(alt, Texo)
        T, number_density = alt_profile[1], alt_profile[8]

        # using eqn(42) of COESA76
        pressure = number_density * k * T
        return pressure

    def density(self, alt, Texo):
        """Solves density at given altitude and exospheric temperature.

        Parameters
        ----------
        alt: ~astropy.units.Quantity
            Geometric/Geopotential altitude.
        Texo: ~astropy.units.Quantity
            Exospheric Temperature

        Returns
        -------
        rho: ~astropy.units.Quantity
            Density at given altitude and exospheric temperature.
        """
        alt_profile = self.altitude_profile(alt, Texo)
        P = self.pressure(alt, Texo)

        # using eqn(42) of COESA76
        rho = P * alt_profile[9] / (R * alt_profile[1])
        return rho
