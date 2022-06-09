import numpy as np
from astropy import units as u
from astropy.coordinates import get_sun
from astropy.time import Time

from poliastro import constants
from poliastro.util import wrap_angle


@u.quantity_input(ltan=u.hourangle)
def raan_from_ltan(epoch, ltan=12.0):
    """RAAN angle from LTAN for SSO around the earth

    Parameters
    ----------
    epoch : ~astropy.time.Time
        Value of time to calculate the RAAN for
    ltan : ~astropy.units.Quantity
        Decimal hour between 0 and 24

    Returns
    -------
    RAAN: ~astropy.units.Quantity
        Right ascension of the ascending node angle in GCRS

    Notes
    -----
    Calculations of the sun mean longitude and equation of time
    follow "Fundamentals of Astrodynamics and Applications"
    Fourth edition by Vallado, David A.
    """

    T_UT1 = ((epoch.ut1 - constants.J2000).value / 36525.0) * u.deg
    T_TDB = ((epoch.tdb - constants.J2000).value / 36525.0) * u.deg

    # Apparent sun position
    sun_position = get_sun(epoch)

    # Calculate the sun apparent local time
    salt = sun_position.ra + 12 * u.hourangle

    # Use the equation of time to calculate the mean sun local time (fictional sun without anomalies)

    # Sun mean anomaly
    M_sun = 357.5291092 * u.deg + 35999.05034 * T_TDB

    # Sun mean longitude
    l_sun = 280.460 * u.deg + 36000.771 * T_UT1
    l_ecliptic_part2 = 1.914666471 * u.deg * np.sin(
        M_sun
    ) + 0.019994643 * u.deg * np.sin(2 * M_sun)
    l_ecliptic = l_sun + l_ecliptic_part2

    eq_time = (
        -l_ecliptic_part2
        + 2.466 * u.deg * np.sin(2 * l_ecliptic)
        - 0.0053 * u.deg * np.sin(4 * l_ecliptic)
    )

    # Calculate sun mean local time

    smlt = salt + eq_time

    # Desired angle between sun and ascending node
    alpha = wrap_angle(ltan, 24 * u.hourangle).to(u.rad)

    # Use the mean sun local time calculate needed RAAN for given LTAN
    raan = smlt + alpha
    return raan


def get_local_sidereal_time(lon, time):
    """Gets the local sideral time from input location.

    Parameters
    ----------
    lon: astropy.units.Quantity
        East longitude of a point on a body.
    time: astropy.time.Time
        Time at which to calculate the local sidereal time.

    Returns
    -------
    theta: astropy.units.Quantity
        Local Sidereal Time (LST).
    """
    current_time = Time(time, scale="utc")

    # Local sidereal time = Greenwich Mean Sidereal time + East longitude of station.
    # theta must be in [0, 360] degree range.
    theta = current_time.sidereal_time("mean", "greenwich") + lon
    if theta < (0 << u.deg) or theta > (360 << u.deg):
        theta = wrap_angle(theta, limit=360 << u.deg)
    return theta
