import numpy as np
from astropy import units as u

pi2 = 1.57079632679
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
R = 8314.32 #* u.J / u.kg / u.K
k = 1.380622e-23
# Na = 6.022169e-26 / u.kmol


class Jacchia77:  # , Z, T, CN2, CO2, CO, CAr, CHe, CH, CM, WM):
    def __init__(self):
        self.E5M = [0.0 for _ in range(11)]
        self.E6P = [0.0 for _ in range(11)]
        self.x = 0.0
        self.y = 0.0
        self.h = 0.0
        self.hbase = 0.0
        self.pbase = 0.0
        self.tbase = 0.0
        self.tgrad = 0.0

    def properties(self, alt, Tinf):
        alt = alt + 1  # in fortran the upper limits are included. in python are not.

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

            # For alt < 86, use U.S. Standard Atmosphere 1976 with added [O].
            if iz < 90:
                print("abhi")
                self.h = self.Z[iz] * 6369.0 / (self.Z[iz] + 6369.0)
                if iz <= 32:
                    if iz <= 11:
                        self.hbase = 0.0
                        self.pbase = 1.0
                        self.tbase = 288.15
                        self.tgrad = -6.5
                        self.goto110(iz)
                        continue
                    elif iz <= 20:
                        self.hbase = 11
                        self.pbase = 2.233611e-1
                        self.tbase = 216.65
                        self.tgrad = 0
                        self.goto120(iz)
                        continue
                    else:
                        self.hbase = 20.0
                        self.pbase = 5.403295e-2
                        self.tbase = 216.65
                        self.tgrad = 1
                        self.goto110(iz)
                        continue
                elif iz <= 51:
                    if iz <= 47:
                        self.hbase = 32.0
                        self.pbase = 8.5666784e-3
                        self.tbase = 228.65
                        self.tgrad = 2.8
                        self.goto110(iz)
                        continue
                    else:
                        self.hbase = 47
                        self.pbase = 1.0945601e-3
                        self.tbase = 270.65
                        self.tgrad = 0
                        self.goto120(iz)
                        continue
                elif iz <= 71:
                    self.hbase = 51.0
                    self.pbase = 6.6063531e-4
                    self.tbase = 270.65
                    self.tgrad = -2.8
                    self.goto110(iz)
                    continue
                else:
                    self.hbase = 71.0
                    self.pbase = 3.9046834e-5
                    self.tbase = 214.65
                    self.tgrad = -2.0
                    self.goto110(iz)
                    continue

            # For 85 < alt < 90, integrate barometric equation with fudged molecular weight
            elif iz <= 89:
                self.T[iz] = 188.0
                self.y = 10.0 ** (
                    -3.7469 + (iz - 85) * (0.226434 - (iz - 85) * 5.945e-3)
                )
                self.WM[iz] = wm0 * (1 - self.y)
                self.CM[iz] = (
                    self.CM[iz - 1]
                    * (self.T[iz - 1] / self.T[iz])
                    * (self.WM[iz] / self.WM[iz - 1])
                    * np.exp(
                        -0.5897446
                        * (
                            (self.WM[iz - 1] / self.T[iz - 1])
                            * (1 + self.Z[iz - 1] / 6356.766) ** (-2)
                            + (self.WM[iz] / self.T[iz])
                            * (1 + self.Z[iz] / 6356.766) ** (-2)
                        )
                    )
                )
                self.goto400(iz)
                continue

            # For alt > 89, use Jacchia 1977
            else:
                if iz <= 90:
                    self.T[iz] = 188
                elif Tinf < 188.1:
                    self.T[iz] = 188
                else:
                    self.x = 0.0045 * (Tinf - 188.0)
                    Tx = 188 + 110.5 * np.log(self.x + np.sqrt(self.x * self.x + 1))
                    Gx = pi2 * 1.9 * (Tx - 188.0) / (125.0 - 90.0)
                    if iz <= 125:
                        self.T[iz] = Tx + ((Tx - 188.0) / pi2) * np.arctan(
                            (Gx / (Tx - 188.0))
                            * (self.Z[iz] - 125.0)
                            * (
                                1.0
                                + 1.7
                                * ((self.Z[iz] - 125.0) / (self.Z[iz] - 90.0)) ** 2
                            )
                        )
                    else:
                        self.T[iz] = Tx + ((Tinf - Tx) / pi2) * np.arctan(
                            (Gx / (Tinf - Tx))
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
                        G0 = (1 + self.Z[iz - 1] / 6356.766) ** (-2)
                        G1 = (1 + self.Z[iz] / 6356.766) ** (-2)
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
                    G0 = (1 + self.Z[iz - 1] / 6356.766) ** (-2)
                    G1 = (1 + self.Z[iz] / 6356.766) ** (-2)

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

                # goto500(alt, Tinf) #These are not needed since they are continues
                continue
            # goto500(alt, Tinf) #These are not needed since they are continues
            # continue

        # return
        # alt_properties = self.goto500(alt, Tinf)
        # return [last for *_, last in alt_properties]
        # if 150 <= alt < 500:
        #     alt_properties = self.goto500(500, Tinf)
        #     dd = [last[alt] for last in alt_properties]
        #     return dd
        #
        #
        #
        # else:
        return self.goto500(alt, Tinf)

    def static_profile(self, alt, Tinf):
        if 150 <= alt < 500:
            alt_properties = self.properties(500, Tinf)
        else:
            alt_properties = self.properties(alt, Tinf)

        return [last[alt] for last in alt_properties]


    def goto110(self, iz):
        self.T[iz] = self.tbase + self.tgrad * (self.h - self.hbase)
        self.x = (self.tbase / self.T[iz]) ** (34.163195 / self.tgrad)
        self.goto130(iz)

    def goto120(self, iz):
        self.T[iz] = self.tbase
        self.x = np.exp(-34.163195 * (self.h - self.hbase) / self.tbase)
        self.goto130(iz)

    def goto130(self, iz):
        self.CM[iz] = 2.547e19 * (288.15 / self.T[iz]) * self.pbase * self.x
        self.goto400(iz)

    def goto400(self, iz):
        # Calculate O/O2 dissociation for Z .lt. 90 km
        self.y = 10.0 ** (-3.7469 + (iz - 85) * (0.226434 - (iz - 85) * 5.945e-3))
        self.x = 1 - self.y
        self.WM[iz] = wm0 * self.x
        self.CN2[iz] = qN2 * self.CM[iz]
        self.CO[iz] = 2.0 * self.y * self.CM[iz]
        self.CO2[iz] = (self.x * qO2 - self.y) * self.CM[iz]
        self.CAr[iz] = qAr * self.CM[iz]
        self.CHe[iz] = qHe * self.CM[iz]
        self.CH[iz] = 0

    def H_correction(self, alt, Tinf):
        phid00 = 10.0 ** (6.9 + 28.9 * Tinf ** (-0.25)) / 2.0e20
        phid00 = phid00 * 5.24e2
        H_500 = 10.0 ** (-0.06 + 28.9 * Tinf ** (-0.25))
        # print(alt)
        for iz in range(150, alt):
            phid0 = phid00 / np.sqrt(self.T[iz])
            self.WM[iz] = (
                    wmH
                    * 0.5897446
                    * ((1.0 + self.Z[iz] / 6356.766) ** (-2))
                    / self.T[iz]
                    + phid0
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
            self.WM[iz] = (wmN2 * self.CN2[iz]
                           + wmO2 * self.CO2[iz]
                           + wmO * self.CO[iz]
                           + wmAr * self.CAr[iz]
                           + wmHe * self.CHe[iz]
                           + wmH * self.CH[iz]
                           ) / self.CM[iz]

    def goto500(self, alt, Tinf):

        # Add Jacchia 1977 empirical corrections to [O] and [O2]
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

        # Calculate [H] from Jacchia 1997 formulas if alt > 500.
        if alt >= 500:
            self.H_correction(alt, Tinf)

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

    def pressure(self, alt, Tinf):
        alt_profile = self.static_profile(alt, Tinf)
        T, number_density = alt_profile[1], alt_profile[8]
        pressure = number_density * k * T
        return pressure#, np.log10(pressure)

    def density(self, alt, Tinf):
        alt_profile = self.static_profile(alt, Tinf)
        P = self.pressure(alt, Tinf)
        rho = P * alt_profile[9] / (R * alt_profile[1])
        return rho, np.log10(rho)







