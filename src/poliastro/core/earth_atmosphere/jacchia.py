"""Low-level calculations for the Jacchia77 atmospheric model.

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
from numba import njit as jit

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
R = 8314.32  # Units: u.J / (u.kg * u.mol)


@jit
def _O_and_O2_correction(alt, Texo, Z, CN2, CO2, CO, CAr, CHe, CH, CM, WM):
    for iz in range(90, alt):
        CO2[iz] = CO2[iz] * (
            10.0 ** (-0.07 * (1.0 + np.tanh(0.18 * (Z[iz] - 111.0))))
        )
        CO[iz] = CO[iz] * (
            10.0 ** (-0.24 * np.exp(-0.009 * (Z[iz] - 97.7) ** 2))
        )
        CM[iz] = CN2[iz] + CO2[iz] + CO[iz] + CAr[iz] + CHe[iz] + CH[iz]
        WM[iz] = (
            wmN2 * CN2[iz]
            + wmO2 * CO2[iz]
            + wmO * CO[iz]
            + wmAr * CAr[iz]
            + wmHe * CHe[iz]
            + wmH * CH[iz]
        ) / CM[iz]


@jit
def _H_correction(alt, Texo, x, y, Z, CN2, CO2, CO, CAr, CHe, CH, CM, WM, T):
    phid00 = 10.0 ** (6.9 + 28.9 * Texo ** (-0.25)) / 2.0e20
    phid00 = phid00 * 5.24e2
    H_500 = 10.0 ** (-0.06 + 28.9 * Texo ** (-0.25))
    # print(alt)
    for iz in range(150, alt):
        phid0 = phid00 / np.sqrt(T[iz])
        WM[iz] = wmH * 0.5897446 * ((1.0 + Z[iz] / R0) ** (-2)) / T[iz] + phid0
        CM[iz] = CM[iz] * phid0

    y = WM[150]
    WM[150] = 0

    for iz in range(151, alt):
        x = WM[iz - 1] + (y + WM[iz])
        y = WM[iz]
        WM[iz] = x

    for iz in range(150, alt):
        WM[iz] = np.exp(WM[iz]) * (T[iz] / T[150]) ** 0.75
        CM[iz] = WM[iz] * CM[iz]

    y = CM[150]
    CM[150] = 0

    for iz in range(151, alt):
        x = CM[iz - 1] + 0.5 * (y + CM[iz])
        y = CM[iz]
        CM[iz] = x

    for iz in range(150, alt):
        CH[iz] = (WM[500] / WM[iz]) * (H_500 - (CM[iz] - CM[500]))

    for iz in range(150, alt):
        CM[iz] = CN2[iz] + CO2[iz] + CO[iz] + CAr[iz] + CHe[iz] + CH[iz]
        WM[iz] = (
            wmN2 * CN2[iz]
            + wmO2 * CO2[iz]
            + wmO * CO[iz]
            + wmAr * CAr[iz]
            + wmHe * CHe[iz]
            + wmH * CH[iz]
        ) / CM[iz]


@jit
def _altitude_profile(alt, Texo, x, y, E5M, E6P):
    # Raise Value Error if alt < 90 km or alt > 2500 km.
    if alt < 90 or 2500 < alt:
        raise ValueError(
            "Jacchia77 has been implemented in range 90km - 2500km."
        )

    alt = int(
        alt + 1
    )  # in fortran the upper limits are included. in python are not.
    Texo = int(Texo)

    Z = [0.0 for _ in range(alt)]
    T = [0.0 for _ in range(alt)]
    CN2 = [0.0 for _ in range(alt)]
    CO2 = [0.0 for _ in range(alt)]
    CO = [0.0 for _ in range(alt)]
    CAr = [0.0 for _ in range(alt)]
    CHe = [0.0 for _ in range(alt)]
    CH = [0.0 for _ in range(alt)]
    CM = [0.0 for _ in range(alt)]
    WM = [0.0 for _ in range(alt)]

    for iz in range(90, alt):
        Z[iz] = iz
        CH[iz] = 0

        if iz <= 90:
            T[iz] = 188
        elif Texo < 188.1:
            T[iz] = 188
        else:
            x = 0.0045 * (Texo - 188.0)
            Tx = 188 + 110.5 * np.log(x + np.sqrt(x * x + 1))
            Gx = pi2 * 1.9 * (Tx - 188.0) / (125.0 - 90.0)
            if iz <= 125:
                T[iz] = Tx + ((Tx - 188.0) / pi2) * np.arctan(
                    (Gx / (Tx - 188.0))
                    * (Z[iz] - 125.0)
                    * (1.0 + 1.7 * ((Z[iz] - 125.0) / (Z[iz] - 90.0)) ** 2)
                )
            else:
                T[iz] = Tx + ((Texo - Tx) / pi2) * np.arctan(
                    (Gx / (Texo - Tx))
                    * (Z[iz] - 125.0)
                    * (1.0 + 5.5e-5 * (Z[iz] - 125.0) ** 2)
                )
        if iz <= 100:
            x = iz - 90
            E5M[iz - 90] = 28.89122 + x * (
                -2.83071e-2
                + x
                * (
                    -6.59924e-3
                    + x * (-3.39574e-4 + x * (+6.19256e-5 + x * (-1.84796e-6)))
                )
            )
            if iz <= 90:
                E6P[0] = 7.145e13 * T[90]
            else:
                G0 = (1 + Z[iz - 1] / R0) ** (-2)
                G1 = (1 + Z[iz] / R0) ** (-2)
                E6P[iz - 90] = E6P[iz - 91] * np.exp(
                    -0.5897446
                    * (
                        G1 * E5M[iz - 90] / T[iz]
                        + G0 * E5M[iz - 91] / T[iz - 1]
                    )
                )

            x = E5M[iz - 90] / wm0
            y = E6P[iz - 90] / T[iz]

            CN2[iz] = qN2 * y * x
            CO[iz] = 2.0 * (1.0 - x) * y
            CO2[iz] = (x * (1.0 + qO2) - 1.0) * y
            CAr[iz] = qAr * y * x
            CHe[iz] = qHe * y * x
            CH[iz] = 0
        else:
            G0 = (1 + Z[iz - 1] / R0) ** (-2)
            G1 = (1 + Z[iz] / R0) ** (-2)

            x = 0.5897446 * (G1 / T[iz] + G0 / T[iz - 1])
            y = T[iz - 1] / T[iz]
            CN2[iz] = CN2[iz - 1] * y * np.exp(-wmN2 * x)
            CO2[iz] = CO2[iz - 1] * y * np.exp(-wmO2 * x)
            CO[iz] = CO[iz - 1] * y * np.exp(-wmO * x)
            CAr[iz] = CAr[iz - 1] * y * np.exp(-wmAr * x)
            CHe[iz] = CHe[iz - 1] * (y**0.62) * np.exp(-wmHe * x)
            CH[iz] = 0

    _O_and_O2_correction(alt, Texo, Z, CN2, CO2, CO, CAr, CHe, CH, CM, WM)

    if 500 <= alt:
        _H_correction(
            alt, Texo, x, y, Z, CN2, CO2, CO, CAr, CHe, CH, CM, WM, T
        )

    return (
        Z,
        T,
        np.array(CN2),
        np.array(CO2),
        np.array(CO),
        np.array(CAr),
        np.array(CHe),
        np.array(CH),
        np.array(CM),
        WM,
    )
