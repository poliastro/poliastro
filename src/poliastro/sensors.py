from astropy import units as u

from poliastro.core.sensors import (
    ground_range_diff_at_azimuth as ground_range_diff_at_azimuth_fast,
    min_and_max_ground_range as min_and_max_ground_range_fast,
)


@u.quantity_input(altitude=u.km, fov=u.rad, boresight=u.rad, R=u.km)
def min_and_max_ground_range(altitude, fov, boresight, R):
    """
    Calculate the minimum and maximum values of ground-range angles.

    Parameters
    ----------
    altitude : ~astropy.units.Quantity
        Altitude over surface.
    fov : ~astropy.units.Quantity
        Angle of the total area that a sensor can observe.
    boresight : ~astropy.units.Quantity
        Center boresight angle.
    R : ~astropy.units.Quantity
        Attractor equatorial radius.

    Returns
    -------
    lat_lon_min: ~astropy.units.Quantity
        Minimum value of latitude and longitude.
    lat_lon_max: ~astropy.units.Quantity
        Maximum value of latitude and longitude.

    """
    altitude = altitude.to_value(u.km)
    fov = fov.to_value(u.rad)
    boresight = boresight.to_value(u.rad)
    R = R.to_value(u.km)
    lat_lon_min, lat_lon_max = min_and_max_ground_range_fast(
        altitude, fov, boresight, R
    )

    return lat_lon_min * u.rad, lat_lon_max * u.rad


@u.quantity_input(
    altitude=u.km,
    fov=u.rad,
    boresight=u.rad,
    azimuth=u.rad,
    nadir_lat=u.rad,
    nadir_lon=u.rad,
    R=u.km,
)
def ground_range_diff_at_azimuth(
    altitude, fov, boresight, azimuth, nadir_lat, nadir_lon, R
):
    """
    Calculate the difference in ground-range angles.

    Use the boresight angle, the latitude and longitude of the target,
    and the desired azimuth (which directs where the sensor is looking).

    Parameters
    ----------
    altitude : ~astropy.units.Quantity
        Altitude over surface.
    fov : ~astropy.units.Quantity
        Angle of the total area that a sensor can observe.
    boresight : ~astropy.units.Quantity
        Center boresight angle.
    azimuth : ~astropy.units.Quantity
        Azimuth angle, used to specify where the sensor is looking.
    nadir_lat : ~astropy.units.Quantity
        Latitude angle of nadir point.
    nadir_lon : ~astropy.units.Quantity
        Longitude angle of nadir point.
    R : ~astropy.units.Quantity
        Attractor equatorial radius.

    Returns
    -------
    ground_range_diff : ~astropy.units.Quantity
        The difference in ground-range angles from the boresight angle.
    target_lat: ~astropy.units.Quantity
        Latitude angle of the target point.
    target_lon: ~astropy.units.Quantity
        Longitude angle of the target point.

    Raises
    ------
    ValueError
        This formula always gives the answer for the short way to the target of the acute azimuth angle,
        which must be greater or equal than 0º and less than 180º.

    """
    altitude = altitude.to_value(u.km)
    fov = fov.to_value(u.rad)
    boresight = boresight.to_value(u.rad)
    azimuth = azimuth.to_value(u.rad)
    nadir_lat = nadir_lat.to_value(u.rad)
    nadir_lon = nadir_lon.to_value(u.rad)
    R = R.to_value(u.km)

    (
        ground_range_diff,
        target_lat,
        target_lon,
    ) = ground_range_diff_at_azimuth_fast(
        altitude, fov, boresight, azimuth, nadir_lat, nadir_lon, R
    )

    return ground_range_diff * u.rad, target_lat * u.rad, target_lon * u.rad
