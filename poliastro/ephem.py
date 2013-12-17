"""Ephemerides.

"""

import numpy as np
from numpy import radians
from numpy.polynomial.polynomial import polyval

from . import angles
from .util import normalize

from astropy import time
from astropy import units
from astropy.constants import au

AU = au.to(units.km).value


__all__ = ['J2000', 'mean_elements']

J2000 = time.Time("2000-01-01 12:00", scale='utc')

MERCURY = 1
VENUS = 2
EARTH = 3
MARS = 4
JUPITER = 5
SATURN = 6
URANUS = 7
NEPTUNE = 8
PLUTO = 9

ephem_coeffs = {
    MERCURY: [
        [0.387098310],
        [0.20563175, 0.000020406, -0.0000000284, -0.00000000017],
        [7.004986, 0.0018215, -0.00001809, 0.000000053],
        [48.330893, 1.1861890, 0.00017587, 0.000000211],
        [77.456119, 1.5564775, 0.00029589, 0.000000056],
        [252.250906, 149474.0722491, 0.00030397, 0.000000018]
    ],
    VENUS: [
        [0.723329820],
        [0.00677188, -0.000047766, 0.0000000975, 0.00000000044],
        [3.394662, 0.0010037, -0.00000088, -0.000000007],
        [76.679920, 0.9011190, 0.00040665, -0.000000080],
        [131.563707, 1.4022188, -0.00107337, -0.000005315],
        [181.979801, 58519.2130302, 0.00031060, 0.000000015]
    ],
    EARTH: [
        [1.000001018],
        [0.01670862, -0.000042037, -0.0000001236, 0.00000000004],
        [0.0],
        [0.0],
        [102.937348, 1.7195269, 0.00045962, 0.000000499],
        [100.4664990, 36000.7698231, 0.00030368, 0.000000021]
    ],
    MARS: [
        [1.523679342],
        [0.09340062, 0.000090483, -0.0000000806, -0.00000000035],
        [1.849726, -0.0006010, 0.00001276, -0.000000006],
        [49.558093, 0.7720923, 0.00001605, 0.000002325],
        [336.060234, 1.8410331, 0.00013515, 0.000000318],
        [355.433275, 19141.6964746, 0.00031097, 0.000000015]
    ],
    JUPITER: [
        [5.202603191, 0.0000001913],
        [0.04849485, 0.000163244, -0.0000004719, -0.00000000197],
        [1.303270, -0.0054966, 0.00000465, -0.000000004],
        [100.464441, 1.0209550, 0.00040117, 0.000000569],
        [14.331309, 1.6126668, 0.00103127, -0.000004569],
        [34.351484, 3036.3027879, 0.00022374, 0.000000025]
    ],
    SATURN: [
        [9.554909596, 0.0000021389],
        [0.05550862, -0.000346818, -0.0000006456, 0.00000000338],
        [2.488878, -0.0037363, -0.00001516, 0.000000089],
        [113.665524, 0.8770970, -0.00012067, -0.000002380],
        [93.056787, 1.9637694, 0.00083757, 0.000004899],
        [50.077471, 1223.5110141, 0.00051952, -0.000000003]
    ],
    URANUS: [
        [19.218446062, -0.0000000372, 0.00000000098],
        [0.04629590, -0.000027337, 0.0000000790, 0.00000000025],
        [0.773196, 0.0007744, 0.00003749, -0.000000092],
        [74.005947, 0.5211258, 0.00133982, 0.000018516],
        [173.005159, 1.4863784, 0.00021450, 0.000000433],
        [314.055005, 429.8640561, 0.00030434, 0.000000026]
    ],
    NEPTUNE: [
        [30.110386869, -0.0000001663,  0.00000000069],
        [0.00898809, 0.000006408, -0.0000000008, -0.00000000005],
        [1.769952, -0.0093082, -0.00000708, 0.000000028],
        [131.784057, 1.1022057, 0.00026006, -0.000000636],
        [48.123691, 1.4262677, 0.00037918, -0.000000003],
        [304.348655, 219.8833092, 0.00030926, 0.000000018]
    ],
    PLUTO: [
        [39.48168677, -0.00076912],
        [0.24880766, 0.00006465],
        [17.13233],
        [110.4065],
        [224.6148],
        [218.88735]
    ]
}


def mean_elements(dd, nbody):
    """Returns the mean orbital elements of a planet for a given date.

    The orbital elements are referred to the mean equator and mean equinox of
    the date.

    Parameters
    ----------
    dd : datetime
        Date to calculate the mean elements.
    nbody : int
        Integer identifying the body.

    Returns
    -------
    a, ecc, inc, omega, argp, nu

    Notes
    -----
    The orbital elements are evaluated as a polynomial, whose oefficients
    are taken from [MEEUS]_.

    References
    ----------
    .. [MEEUS] Meeus, "Astronomical Algorithms", Willmann-Bell, Incorporated,
       pp. 200-202, 1991.

    Examples
    --------
    >>> from datetime import datetime
    >>> from poliastro import ephem
    >>> dd = datetime(2065, 6, 24)
    >>> ephem.mean_elements(dd, ephem.MERCURY)
    array([  5.79090829e+07,   2.05645099e-01,   1.22280751e-01,
             8.57090185e-01,   5.12563608e-01,   2.47136342e+00])

    TODO: Check with ASTINTER.for.

    """
    coeffs = ephem_coeffs[nbody]
    tt = (time.Time(dd, scale='utc') - J2000).jd / 36525
    a = polyval(tt, coeffs[0]) * AU
    ecc = polyval(tt, coeffs[1])
    inc = normalize(radians(polyval(tt, coeffs[2])))
    omega = radians(polyval(tt, coeffs[3]))
    lonper = radians(polyval(tt, coeffs[4]))
    meanlon = radians(polyval(tt, coeffs[5]))
    argp = normalize(lonper - omega)
    M = normalize(meanlon - lonper)
    _, nu = angles.M2nu(ecc, M)
    coe = np.column_stack((a, ecc, inc, normalize(omega), argp, nu))
    if coe.shape[0] == 1:
        coe = coe[0]
    return coe
