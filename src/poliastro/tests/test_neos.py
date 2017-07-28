from unittest import mock
import requests

import astropy.units as u
from poliastro import neos

@mock.patch('poliastro.neos.requests.get')
def test_orbit_from_spk_id_has_proper_values(mock_get):

    mock_response = mock.Mock(requests.Response())

    mock_orbital_data = {
        'orbital_data': {
            'orbit_id': '611',
            'orbit_determination_date': '2017-06-06 06:20:43',
            'orbit_uncertainty': '0',
            'minimum_orbit_intersection': '.150505',
            'jupiter_tisserand_invariant': '4.583',
            'epoch_osculation': '2458000.5',
            'eccentricity': '.2225889698301071',
            'semi_major_axis': '1.457940027185708',
            'inclination': '10.82759100494802',
            'ascending_node_longitude': '304.3221633898424',
            'orbital_period': '642.9954742523677',
            'perihelion_distance': '1.133418658460363',
            'perihelion_argument': '178.8165910886752',
            'aphelion_distance': '1.782461395911054',
            'perihelion_time': '2457873.186399333365',
            'mean_anomaly': '71.28027812836476',
            'mean_motion': '.5598795239089109',
            'equinox': 'J2000'
        }
    }

    # mock_response.status_code = 200
    mock_response.json.return_value = mock_orbital_data
    mock_get.return_value = mock_response
    ss = neos.get_orbit_from_spk_id('')

    assert ss.a == mock_orbital_data['orbital_data']['semi_major_axis'] * u.AU
