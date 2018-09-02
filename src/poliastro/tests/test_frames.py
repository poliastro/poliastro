import pytest

from astropy import units as u
from astropy.coordinates import (
    CartesianRepresentation,
    get_body_barycentric, solar_system_ephemeris
)
from astropy.tests.helper import assert_quantity_allclose

from poliastro.constants import J2000
from poliastro.bodies import (
    Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto,
)
from poliastro.frames import (
    ICRS,
    HCRS, MercuryICRS, VenusICRS, GCRS, MarsICRS, JupiterICRS, SaturnICRS, UranusICRS, NeptuneICRS, PlutoICRS
)


@pytest.mark.parametrize("body, frame", [
    (Sun, HCRS),
    (Mercury, MercuryICRS),
    (Venus, VenusICRS),
    (Earth, GCRS),
    (Mars, MarsICRS),
    (Jupiter, JupiterICRS),
    (Saturn, SaturnICRS),
    (Uranus, UranusICRS),
    (Neptune, NeptuneICRS),
    (Pluto, PlutoICRS),
])
def test_planetary_icrs_frame_is_just_translation(body, frame):
    with solar_system_ephemeris.set("de432s"):
        epoch = J2000
        vector = CartesianRepresentation(x=100 * u.km, y=100 * u.km, z=100 * u.km)
        vector_result = frame(vector, obstime=epoch).transform_to(ICRS).represent_as(CartesianRepresentation)

        expected_result = get_body_barycentric(body.name, epoch) + vector

    assert_quantity_allclose(vector_result.xyz, expected_result.xyz)


@pytest.mark.parametrize("body, frame", [
    (Sun, HCRS),
    (Mercury, MercuryICRS),
    (Venus, VenusICRS),
    (Earth, GCRS),
    (Mars, MarsICRS),
    (Jupiter, JupiterICRS),
    (Saturn, SaturnICRS),
    (Uranus, UranusICRS),
    (Neptune, NeptuneICRS),
    (Pluto, PlutoICRS),
])
def test_icrs_body_position_to_planetary_frame_yields_zeros(body, frame):
    with solar_system_ephemeris.set("de432s"):
        epoch = J2000
        vector = get_body_barycentric(body.name, epoch)

        vector_result = ICRS(vector).transform_to(frame(obstime=epoch)).represent_as(CartesianRepresentation)

    assert_quantity_allclose(vector_result.xyz, [0, 0, 0] * u.km, atol=1e-7 * u.km)
