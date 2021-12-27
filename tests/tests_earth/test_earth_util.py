from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time

from poliastro.earth.util import raan_from_ltan


def test_raan_from_ltan_metopb():
    # MetOp-B LTAN: 21:31:45
    # LTAN from https://www.ospo.noaa.gov/Operations/METOP/status.html
    # METOP-B
    # 1 38771U 12049A   20049.95408566 -.00000014  00000-0  13607-4 0  9997
    # 2 38771  98.7092 110.9899 0001263  48.5458 295.8781 14.21485930385043

    ltan = (21 + ((31 + 45 / 60) / 60)) * u.hourangle
    epoch = Time(
        Time("2020-01-01 00:00").to_value("mjd") + 49.95408566 - 1,
        format="mjd",
    )
    expected_raan = 110.9899 * u.deg

    raan = raan_from_ltan(epoch, ltan)

    assert_quantity_allclose(raan.wrap_at(360 * u.deg), expected_raan, atol=0.3 * u.deg)


def test_raan_from_ltan_sentinel5p():
    # SENTINEL-5P LTAN: 13:30
    # LTAN from https://sentinels.copernicus.eu/web/sentinel/missions/sentinel-5p/geographical-coverage
    # 1 42969U 17064A   20049.78099017 -.00000032  00000-0  54746-5 0  9991
    # 2 42969  98.7249 350.5997 0001077  82.0109 278.1189 14.19549365121775

    ltan = (13 + (30 / 60)) * u.hourangle
    epoch = Time(
        Time("2020-01-01 00:00").to_value("mjd") + 49.78099017 - 1,
        format="mjd",
    )
    expected_raan = 350.5997 * u.deg

    raan = raan_from_ltan(epoch, ltan)

    assert_quantity_allclose(raan.wrap_at(360 * u.deg), expected_raan, atol=0.3 * u.deg)
