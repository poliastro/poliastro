import os.path
from unittest import mock

import astropy.units as u
import pytest
import requests
from astropy.tests.helper import assert_quantity_allclose

from poliastro.frames import HeliocentricEclipticJ2000
from poliastro.neos import neows
from poliastro.twobody.angles import nu_to_M

CURRENT_DIRECTORY = os.path.dirname(os.path.realpath(__file__))


@mock.patch("poliastro.neos.neows.requests.Response")
@mock.patch("poliastro.neos.neows.requests.get")
def test_orbit_from_spk_id_has_proper_values(mock_get, mock_response):
    mock_orbital_data = {
        "orbital_data": {
            "eccentricity": ".2225889698301071",
            "semi_major_axis": "1.457940027185708",
            "inclination": "10.82759100494802",
            "ascending_node_longitude": "304.3221633898424",
            "perihelion_argument": "178.8165910886752",
            "mean_anomaly": "71.28027812836476",
            "epoch_osculation": "2458000.5",
        }
    }

    mock_response.json.return_value = mock_orbital_data
    mock_get.return_value = mock_response
    ss = neows.orbit_from_spk_id("")

    assert ss.frame.is_equivalent_frame(HeliocentricEclipticJ2000(obstime=ss.epoch))
    assert ss.ecc == float(mock_orbital_data["orbital_data"]["eccentricity"]) * u.one
    assert ss.a == float(mock_orbital_data["orbital_data"]["semi_major_axis"]) * u.AU
    assert ss.inc == float(mock_orbital_data["orbital_data"]["inclination"]) * u.deg
    assert (
        ss.raan
        == float(mock_orbital_data["orbital_data"]["ascending_node_longitude"]) * u.deg
    )
    assert (
        ss.argp
        == float(mock_orbital_data["orbital_data"]["perihelion_argument"]) * u.deg
    )
    assert_quantity_allclose(
        nu_to_M(ss.nu, ss.ecc),
        float(mock_orbital_data["orbital_data"]["mean_anomaly"]) * u.deg,
        rtol=1e-8,
    )


@mock.patch("poliastro.neos.neows.requests.get")
def test_orbit_from_spk_id_raises_when_error(mock_get):
    resp = requests.Response()

    resp.status_code = 404
    mock_get.return_value = resp
    with pytest.raises(requests.HTTPError):
        neows.orbit_from_spk_id("")


@mock.patch("poliastro.neos.neows.requests.get")
def test_spk_id_from_name_raises_when_error(mock_get):
    resp = requests.Response()

    resp.status_code = 404
    mock_get.return_value = resp
    with pytest.raises(requests.HTTPError):
        neows.spk_id_from_name("")


@mock.patch("poliastro.neos.neows.requests.Response")
@mock.patch("poliastro.neos.neows.requests.get")
def test_spk_id_from_name_parses_body(mock_get, mock_response):
    with open(os.path.join(CURRENT_DIRECTORY, "table.html"), "r") as demo_html:
        html = demo_html.read().replace("\n", "")

    mock_response.text = html
    mock_get.return_value = mock_response
    assert "2000433" == neows.spk_id_from_name("")


@mock.patch("poliastro.neos.neows.requests.Response")
@mock.patch("poliastro.neos.neows.requests.get")
def test_spk_id_from_name_parses_object_list_and_raises(mock_get, mock_response):
    with open(os.path.join(CURRENT_DIRECTORY, "center.html"), "r") as demo_html:
        html = demo_html.read().replace("\n", "")

    mock_response.text = html
    mock_get.return_value = mock_response
    with pytest.raises(ValueError) as e_msg:
        neows.spk_id_from_name("")
        assert "different bodies found" in str(e_msg)


@mock.patch("poliastro.neos.neows.requests.Response")
@mock.patch("poliastro.neos.neows.requests.get")
def test_spk_id_from_name_raises_when_not_found(mock_get, mock_response):
    with open(os.path.join(CURRENT_DIRECTORY, "none.html"), "r") as demo_html:
        html = demo_html.read().replace("\n", "")
    mock_response.text = html
    mock_get.return_value = mock_response
    with pytest.raises(ValueError) as e_msg:
        neows.spk_id_from_name("")
        assert "Object could not be found" in str(e_msg)
