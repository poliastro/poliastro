"""
Low level ra, dec, and W calculations for bodies.

REFERENCES:
The values for all bodies except the Moon have been taken from: Archinal, B., \
Acton, C., A’Hearn, M., Conrad, A., Consolmagno, G., & Duxbury, T. et al. \
(2018). Report of the IAU Working Group on Cartographic Coordinates and \
Rotational Elements: 2015. Celestial Mechanics And Dynamical Astronomy, 130(3). \
doi: 10.1007/s10569-017-9805-5.

The values for Moon have been taken from: Archinal, B., A’Hearn, M., Bowell, E., \
Conrad, A., Consolmagno, G., & Courtin, R. et al. (2010). Report of the IAU \
Working Group on Cartographic Coordinates and Rotational Elements: 2009. \
Celestial Mechanics And Dynamical Astronomy, 109(2), 101-135. doi: 10.1007/s10569-010-9320-4.

"""

import numpy as np
from numba import njit as jit


@jit
def sun_rot_elements_at_epoch(T, d):
    """Calculate rotational elements for Sun.

    Parameters
    ----------
    T : float
        Interval from the standard epoch, in Julian centuries i.e. 36525 days.
    d : float
        Interval in days from the standard epoch.

    Returns
    -------
    ra, dec, W: tuple (float)
        Right ascension and declination of north pole, and angle of the prime meridian.

    """
    ra = 286.13
    dec = 63.87
    W = 84.176 + 14.1844000 * d

    return ra, dec, W


@jit
def mercury_rot_elements_at_epoch(T, d):
    """Calculate rotational elements for Mercury.

    Parameters
    ----------
    T : float
        Interval from the standard epoch, in Julian centuries.
    d : float
        Interval in days from the standard epoch.

    Returns
    -------
    ra, dec, W: tuple (float)
        Right ascension and declination of north pole, and angle of the prime meridian.

    """
    M1 = np.deg2rad(174.7910857 + 4.092335 * d)
    M2 = np.deg2rad(349.5821714 + 8.184670 * d)
    M3 = np.deg2rad(164.3732571 + 12.277005 * d)
    M4 = np.deg2rad(339.1643429 + 16.369340 * d)
    M5 = np.deg2rad(153.9554286 + 20.461675 * d)
    ra = 281.0103 - 0.0328 * T
    dec = 61.45 - 0.005 * T
    W = (329.5988 + 6.1385108 * d) + (
        0.01067257 * np.sin(M1)
        - 0.00112309 * np.sin(M2)
        - 0.00011040 * np.sin(M3)
        - 0.00002539 * np.sin(M4)
        - 0.00000571 * np.sin(M5)
    )

    return ra, dec, W


@jit
def venus_rot_elements_at_epoch(T, d):
    """Calculate rotational elements for Venus.

    Parameters
    ----------
    T : float
        Interval from the standard epoch, in Julian centuries.
    d : float
        Interval in days from the standard epoch.

    Returns
    -------
    ra, dec, W: tuple (float)
        Right ascension and declination of north pole, and angle of the prime meridian.

    """
    ra = 272.76
    dec = 67.16
    W = 160.20 - 1.4813688 * d

    return ra, dec, W


@jit
def mars_rot_elements_at_epoch(T, d):
    """Calculate rotational elements for Mars.

    Parameters
    ----------
    T : float
        Interval from the standard epoch, in Julian centuries
    d : float
        Interval in days from the standard epoch

    Returns
    -------
    ra, dec, W: tuple (float)
        Right ascension and declination of north pole, and angle of the prime meridian.

    """
    M1 = np.deg2rad(198.991226 + 19139.4819985 * T)
    M2 = np.deg2rad(226.292679 + 38280.8511281 * T)
    M3 = np.deg2rad(249.663391 + 57420.7251593 * T)
    M4 = np.deg2rad(266.183510 + 76560.6367950 * T)
    M5 = np.deg2rad(79.398797 + 0.5042615 * T)

    ra = (
        317.269202
        - 0.10927547 * T
        + 0.000068 * np.sin(M1)
        + 0.000238 * np.sin(M2)
        + 0.000052 * np.sin(M3)
        + 0.000009 * np.sin(M4)
        + 0.419057 * np.sin(M5)
    )

    K1 = np.deg2rad(122.433576 + 19139.9407476 * T)
    K2 = np.deg2rad(43.058401 + 38280.8753272 * T)
    K3 = np.deg2rad(57.663379 + 57420.7517205 * T)
    K4 = np.deg2rad(79.476401 + 76560.6495004 * T)
    K5 = np.deg2rad(166.325722 + 0.5042615 * T)

    dec = (
        54.432516
        - 0.05827105 * T
        + 0.000051 * np.cos(K1)
        + 0.000141 * np.cos(K2)
        + 0.000031 * np.cos(K3)
        + 0.000005 * np.cos(K4)
        + 1.591274 * np.cos(K5)
    )

    J1 = np.deg2rad(129.071773 + 19140.0328244 * T)
    J2 = np.deg2rad(36.352167 + 38281.0473591 * T)
    J3 = np.deg2rad(56.668646 + 57420.9295360 * T)
    J4 = np.deg2rad(67.364003 + 76560.2552215 * T)
    J5 = np.deg2rad(104.792680 + 95700.4387578 * T)
    J6 = np.deg2rad(95.391654 + 0.5042615 * T)

    W = (
        176.049863
        + 350.891982443297 * d
        + 0.000145 * np.sin(J1)
        + 0.000157 * np.sin(J2)
        + 0.000040 * np.sin(J3)
        + 0.000001 * np.sin(J4)
        + 0.000001 * np.sin(J5)
        + 0.584542 * np.sin(J6)
    )

    return ra, dec, W


@jit
def jupiter_rot_elements_at_epoch(T, d):
    """Calculate rotational elements for Jupiter.

    Parameters
    ----------
    T : float
        Interval from the standard epoch, in Julian centuries
    d : float
        Interval in days from the standard epoch

    Returns
    -------
    ra, dec, W: tuple (float)
        Right ascension and declination of north pole, and angle of the prime meridian.

    """
    Ja = np.deg2rad(99.360714 + 4850.4046 * T)
    Jb = np.deg2rad(175.895369 + 1191.9605 * T)
    Jc = np.deg2rad(300.323162 + 262.5475 * T)
    Jd = np.deg2rad(114.012305 + 6070.2476 * T)
    Je = np.deg2rad(49.511251 + 64.3000 * T)

    ra = (
        268.056595
        - 0.006499 * T
        + 0.000117 * np.sin(Ja)
        + 0.000938 * np.sin(Jb)
        + 0.001432 * np.sin(Jc)
        + 0.000030 * np.sin(Jd)
        + 0.002150 * np.sin(Je)
    )
    dec = (
        64.495303
        + 0.002413 * T
        + 0.000050 * np.cos(Ja)
        + 0.000404 * np.cos(Jb)
        + 0.000617 * np.cos(Jc)
        - 0.000013 * np.cos(Jd)
        + 0.000926 * np.cos(Je)
    )
    W = 284.95 + 870.536 * d

    return ra, dec, W


@jit
def saturn_rot_elements_at_epoch(T, d):
    """Calculate rotational elements for Saturn.

    Parameters
    ----------
    T : float
        Interval from the standard epoch, in Julian centuries.
    d : float
        Interval in days from the standard epoch.

    Returns
    -------
    ra, dec, W: tuple (float)
        Right ascension and declination of north pole, and angle of the prime meridian.

    """
    ra = 40.589 - 0.036 * T
    dec = 83.537 - 0.004 * T
    W = 38.90 + 810.7939024 * d

    return ra, dec, W


@jit
def uranus_rot_elements_at_epoch(T, d):
    """Calculate rotational elements for Uranus.

    Parameters
    ----------
    T : float
        Interval from the standard epoch, in Julian centuries.
    d : float
        Interval in days from the standard epoch.

    Returns
    -------
    ra, dec, W: tuple (float)
        Right ascension and declination of north pole, and angle of the prime meridian.

    """
    ra = 257.311
    dec = -15.175
    W = 203.81 - 501.1600928 * d

    return ra, dec, W


@jit
def neptune_rot_elements_at_epoch(T, d):
    """Calculate rotational elements for Neptune.

    Parameters
    ----------
    T : float
        Interval from the standard epoch, in Julian centuries.
    d : float
        Interval in days from the standard epoch.

    Returns
    -------
    ra, dec, W: tuple (float)
        Right ascension and declination of north pole, and angle of the prime meridian.

    """
    N = np.deg2rad(357.85 + 52.316 * T)

    ra = 299.36 + 0.70 * np.sin(N)
    dec = 43.46 - 0.51 * np.cos(N)
    W = 249.978 + 541.1397757 * d - 0.48 * np.sin(N)

    return ra, dec, W


@jit
def moon_rot_elements_at_epoch(T, d):
    """Calculate rotational elements for Moon.

    Parameters
    ----------
    T : float
        Interval from the standard epoch, in Julian centuries.
    d : float
        Interval in days from the standard epoch.

    Returns
    -------
    ra, dec, W: tuple (float)
        Right ascension and declination of north pole, and angle of the prime meridian.

    """
    E1 = np.deg2rad(125.045 - 0.0529921 * d)
    E2 = np.deg2rad(250.089 - 0.1059842 * d)
    E3 = np.deg2rad(260.008 + 13.0120009 * d)
    E4 = np.deg2rad(176.625 + 13.3407154 * d)
    E5 = np.deg2rad(357.529 + 0.9856003 * d)
    E6 = np.deg2rad(311.589 + 26.4057084 * d)
    E7 = np.deg2rad(134.963 + 13.0649930 * d)
    E8 = np.deg2rad(276.617 + 0.3287146 * d)
    E9 = np.deg2rad(34.226 + 1.7484877 * d)
    E10 = np.deg2rad(15.134 - 0.1589763 * d)
    E11 = np.deg2rad(119.743 + 0.0036096 * d)
    E12 = np.deg2rad(239.961 + 0.1643573 * d)
    E13 = np.deg2rad(25.053 + 12.9590088 * d)

    ra = (
        269.9949
        + 0.0031 * T
        - 3.8787 * np.sin(E1)
        - 0.1204 * np.sin(E2)
        + 0.0700 * np.sin(E3)
        - 0.0172 * np.sin(E4)
        + 0.0072 * np.sin(E6)
        - 0.0052 * np.sin(E10)
        + 0.0043 * np.sin(E13)
    )
    dec = (
        66.5392
        + 0.0130 * T
        + 1.5419 * np.cos(E1)
        + 0.0239 * np.cos(E2)
        - 0.0278 * np.cos(E3)
        + 0.0068 * np.cos(E4)
        - 0.0029 * np.cos(E6)
        + 0.0009 * np.cos(E7)
        + 0.0008 * np.cos(E10)
        - 0.0009 * np.cos(E13)
    )
    W = (
        38.3213
        + 13.17635815 * d
        - 1.4e-12 * d**2
        + 3.5610 * np.sin(E1)
        + 0.1208 * np.sin(E2)
        - 0.0642 * np.sin(E3)
        + 0.0158 * np.sin(E4)
        + 0.0252 * np.sin(E5)
        - 0.0066 * np.sin(E6)
        - 0.0047 * np.sin(E7)
        - 0.0046 * np.sin(E8)
        + 0.0028 * np.sin(E9)
        + 0.0052 * np.sin(E10)
        + 0.0040 * np.sin(E11)
        + 0.0019 * np.sin(E12)
        - 0.0044 * np.sin(E13)
    )

    return ra, dec, W
