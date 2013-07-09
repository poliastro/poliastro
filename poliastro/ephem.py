"""Ephemerides.

Function gcal2jd borrowed from jdcal by Prasanth Nair, available under
the same license as poliastro (BSD). See COPYING for full text.

"""

import numpy as np
from numpy import radians
from numpy.polynomial.polynomial import polyval

from . import angles
from .constants import AU


__all__ = ['J2000', 'mean_elements']

MJD_0 = 2400000.5
MJD_JD2000 = 51544.5
J2000 = MJD_0 + MJD_JD2000

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
    dd : date
        Date of the query.
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
    >>> ephem.ephem(dd, ephem.MERCURY)
    (57909082.92756851, 0.20564509902750358, 0.12228075058282453, 0.85709018507404022, 0.51256360806673706, 2.4713634244634277)

    """
    coeffs = ephem_coeffs[nbody]
    tt = (jd(dd) - J2000) / 36525
    a = polyval(tt, coeffs[0]) * AU
    ecc = polyval(tt, coeffs[1])
    inc = radians(polyval(tt, coeffs[2]))
    omega = radians(polyval(tt, coeffs[3]))
    lonper = radians(polyval(tt, coeffs[4]))
    meanlon = radians(polyval(tt, coeffs[5]))
    argp = lonper - omega
    M = meanlon - lonper
    _, nu = angles.M2nu(ecc, M)
    return a, ecc, inc, omega, argp, nu


def jd(dd):
    """Returns Julian Day given date.

    Parameters
    ----------
    dd : date or datetime
        Date to compute the JD (Gregorian calendar)

    Examples
    --------
    >>> from datetime import date, datetime
    >>> from poliastro import ephem
    >>> dd = date(2000, 1, 1)
    >>> ephem.jd(dd)
    2451545.0
    >>> ephem.jd(datetime(2065, 6, 24, 0))
    2475460.5
    >>> ephem.jd(datetime(1977, 4, 26, 9, 36))
    2443259.9
    >>> ephem.jd(datetime(1957, 10, 4, 19, 26, 24))
    2436116.31

    """
    try:
        tt = dd.timetz()
        if tt.tzinfo:
            raise NotImplementedError("Only UTC accepted")
        frac = (3600000000 * tt.hour +
                60000000 * tt.minute +
                1000000 * tt.second +
                tt.microsecond) / (24 * 3600 * 1000000)
    except AttributeError:
        # The object is a date: assume midnight
        frac = 0.5
    return sum(gcal2jd(dd.year, dd.month, dd.day) + (frac,))


def gcal2jd(year, month, day):
    """Gregorian calendar date to Julian date.

    The input and output are for the proleptic Gregorian calendar,
    i.e., no consideration of historical usage of the calendar is
    made.

    # TODO: Vectorize, merge with jd?

    Parameters
    ----------
    year : int
        Year as an integer.
    month : int
        Month as an integer.
    day : int
        Day as an integer.

    Returns
    -------
    jd1, jd2: 2-element tuple of floats
        When added together, the numbers give the Julian date for the
        given Gregorian calendar date. The first number is always
        MJD_0 i.e., 2451545.5. So the second is the MJD.

    Examples
    --------
    >>> gcal2jd(2000,1,1)
    (2400000.5, 51544.0)
    >>> 2400000.5 + 51544.0 + 0.5
    2451545.0
    >>> year = [-4699, -2114, -1050, -123, -1, 0, 1, 123, 1678.0, 2000,
    ....: 2012, 2245]
    >>> month = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    >>> day = [1, 12, 23, 14, 25, 16, 27, 8, 9, 10, 11, 31]
    >>> x = [gcal2jd(y, m, d) for y, m, d in zip(year, month, day)]
    >>> for i in x: print(i)
    (2400000.5, -2395215.0)
    (2400000.5, -1451021.0)
    (2400000.5, -1062364.0)
    (2400000.5, -723762.0)
    (2400000.5, -679162.0)
    (2400000.5, -678774.0)
    (2400000.5, -678368.0)
    (2400000.5, -633797.0)
    (2400000.5, -65812.0)
    (2400000.5, 51827.0)
    (2400000.5, 56242.0)
    (2400000.5, 141393.0)

    Negative months and days are valid. For example, 2000/-2/-4 =>
    1999/+12-2/-4 => 1999/10/-4 => 1999/9/30-4 => 1999/9/26.

    >>> gcal2jd(2000, -2, -4)
    (2400000.5, 51447.0)
    >>> gcal2jd(1999, 9, 26)
    (2400000.5, 51447.0)

    >>> gcal2jd(2000, 2, -1)
    (2400000.5, 51573.0)
    >>> gcal2jd(2000, 1, 30)
    (2400000.5, 51573.0)

    >>> gcal2jd(2000, 3, -1)
    (2400000.5, 51602.0)
    >>> gcal2jd(2000, 2, 28)
    (2400000.5, 51602.0)

    Month 0 becomes previous month.

    >>> gcal2jd(2000, 0, 1)
    (2400000.5, 51513.0)
    >>> gcal2jd(1999, 12, 1)
    (2400000.5, 51513.0)

    Day number 0 becomes last day of previous month.

    >>> gcal2jd(2000, 3, 0)
    (2400000.5, 51603.0)
    >>> gcal2jd(2000, 2, 29)
    (2400000.5, 51603.0)

    If `day` is greater than the number of days in `month`, then it
    gets carried over to the next month.

    >>> gcal2jd(2000,2,30)
    (2400000.5, 51604.0)
    >>> gcal2jd(2000,3,1)
    (2400000.5, 51604.0)

    >>> gcal2jd(2001,2,30)
    (2400000.5, 51970.0)
    >>> gcal2jd(2001,3,2)
    (2400000.5, 51970.0)

    Notes
    -----
    The returned Julian date is for mid-night of the given date. To
    find the Julian date for any time of the day, simply add time as a
    fraction of a day. For example Julian date for mid-day can be
    obtained by adding 0.5 to either the first part or the second
    part. The latter is preferable, since it will give the MJD for the
    date and time.

    BC dates should be given as -(BC - 1) where BC is the year. For
    example 1 BC == 0, 2 BC == -1, and so on.

    Negative numbers can be used for `month` and `day`. For example
    2000, -1, 1 is the same as 1999, 11, 1.

    The Julian dates are proleptic Julian dates, i.e., values are
    returned without considering if Gregorian dates are valid for the
    given date.

    The input values are truncated to integers.

    """
    year = int(year)
    month = int(month)
    day = int(day)

    a = np.trunc((month - 14) / 12.0)
    jd = np.trunc((1461 * (year + 4800 + a)) / 4.0)
    jd += np.trunc((367 * (month - 2 - 12 * a)) / 12.0)
    x = np.trunc((year + 4900 + a) / 100.0)
    jd -= np.trunc((3 * x) / 4.0)
    jd += day - 2432075.5  # was 32075; add 2400000.5

    jd -= 0.5  # 0 hours; above JD is for midday, switch to midnight.

    return MJD_0, jd
