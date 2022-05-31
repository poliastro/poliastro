import numpy as np
from astropy import units as u

from poliastro.core.earth_atmosphere.jacchia import (
    _altitude_profile as _altitude_profile_fast,
    _H_correction as _H_correction_fast,
    _O_and_O2_correction as _O_and_O2_correction_fast,
    wmAr,
    wmH,
    wmHe,
    wmN2,
    wmO,
    wmO2,
)

R = 8314.32 * u.J / (u.K * u.kmol)
k = 1.380622e-23 * u.J / u.K
Na = R / k


class Jacchia77:
    """Holds the model for U.S Standard Atmosphere 1962."""

    def __init__(self, Texo):
        self.E5M = np.zeros(11)
        self.E6P = np.zeros(11)
        self.x = 0.0
        self.y = 0.0
        self.Texo = Texo

    def _altitude_profile(self, alt):
        Z, T, CN2, CO2, CO, CAr, CHe, CH, CM, WM = _altitude_profile_fast(
            alt, self.Texo.to_value(u.K), self.x, self.y, self.E5M, self.E6P
        )
        (
            self.Z,
            self.T,
            self.CN2,
            self.CO2,
            self.CO,
            self.CAr,
            self.CHe,
            self.CH,
            self.CM,
            self.WM,
        ) = (
            Z * u.km,
            T * u.K,
            np.array(CN2) * 1e6 * (u.m) ** -3,
            np.array(CO2) * 1e6 * (u.m) ** -3,
            np.array(CO) * 1e6 * (u.m) ** -3,
            np.array(CAr) * 1e6 * (u.m) ** -3,
            np.array(CHe) * 1e6 * (u.m) ** -3,
            np.array(CH) * 1e6 * (u.m) ** -3,
            np.array(CM) * 1e6 * (u.m) ** -3,
            np.array(WM) * 1e-3 * (u.kg / u.mol),
        )
        return (
            self.Z,
            self.T,
            self.CN2,
            self.CO2,
            self.CO,
            self.CAr,
            self.CHe,
            self.CH,
            self.CM,
            self.WM,
        )

    def _H_correction(self, alt):
        """Calculate [H] from Jacchia 1977 formulas"""
        _H_correction_fast(alt, self.Texo)

    def _O_and_O2_correction(self, alt):
        """Add Jacchia 1977 empirical corrections to [O] and [O2]"""
        _O_and_O2_correction_fast(alt, self.Texo)

    def altitude_profile(self, alt):
        """Solves for atmospheric altitude profile at given altitude and exospheric temperature.

        Parameters
        ----------
        alt : ~astropy.units.Quantity
            Geometric/Geopotential altitude.

        Returns
        -------
        altitude_profile: list
            [altitude(Z), T, N2, O2, O, Ar, He, H, Total number density, Mean Molecular weight]
        """
        # checking if the units entered are km
        if alt.unit == u.km:
            if 150 <= alt.value < 500:
                alt_properties = self._altitude_profile(500)
            else:
                alt_properties = self._altitude_profile(alt.value)

        return [last[int(alt.value)] for last in alt_properties]

    def temperature(self, alt):
        """Solves for temperature at given altitude and exospheric temperature.

        Parameters
        ----------
        alt : ~astropy.units.Quantity
            Geometric/Geopotential altitude.

        Returns
        -------
        T: ~astropy.units.Quantity
            Absolute temeperature  and exospheric temperature
        """
        T = self.altitude_profile(alt)[1]
        return T

    def pressure(self, alt):
        """Solves pressure at given altitude and exospheric temperature.

        Parameters
        ----------
        alt : ~astropy.units.Quantity
            Geometric/Geopotential altitude.

        Returns
        -------
        p: ~astropy.units.Quantity
            Pressure at given altitude  and exospheric temperature.
        """
        alt_profile = self.altitude_profile(alt)
        T, number_density = alt_profile[1], alt_profile[8]

        # using eqn(42) of COESA76
        pressure = number_density * k * T
        return pressure

    def density(self, alt):
        """Solves density at given altitude and exospheric temperature.

        Parameters
        ----------
        alt : ~astropy.units.Quantity
            Geometric/Geopotential altitude.

        Returns
        -------
        rho: ~astropy.units.Quantity
            Density at given altitude and exospheric temperature.
        """
        (Z, T, CN2, CO2, CO, CAr, CHe, CH, CM, WM) = self.altitude_profile(alt)

        # using eqn(42) of COESA for multiple gases
        M_i = [wmN2, wmO2, wmO, wmAr, wmHe, wmH] << (u.g / u.mol)
        n_i = [
            CN2.to_value(u.m**-3),
            CO2.to_value(u.m**-3),
            CO.to_value(u.m**-3),
            CAr.to_value(u.m**-3),
            CHe.to_value(u.m**-3),
            CH.to_value(u.m**-3),
        ] << (1 / u.m**3)
        rho = (n_i @ M_i) / Na
        return rho.to(u.kg / u.m**3)
