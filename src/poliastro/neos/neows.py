"""NEOs orbit from NEOWS and JPL SBDB

"""
import re

from bs4 import BeautifulSoup
import requests

import astropy.units as u
from astropy.time import Time

from poliastro.twobody.orbit import Orbit
from poliastro.bodies import Sun
from poliastro.frames import Planes
from poliastro.twobody.angles import M_to_nu

# Base URLs
NEOWS_URL = 'https://api.nasa.gov/neo/rest/v1/neo/'
SBDB_URL = 'https://ssd.jpl.nasa.gov/sbdb.cgi'
DEFAULT_API_KEY = 'DEMO_KEY'


def orbit_from_spk_id(spk_id, api_key=None):
    """Return :py:class:`~poliastro.twobody.orbit.Orbit` given a SPK-ID.

    Retrieve info from NASA NeoWS API, and therefore
    it only works with NEAs (Near Earth Asteroids).

    Parameters
    ----------
    spk_id : str
        SPK-ID number, which is given to each body by JPL.
    api_key : str
        NASA OPEN APIs key (default: `DEMO_KEY`)

    Returns
    -------
    orbit : ~poliastro.twobody.orbit.Orbit
        NEA orbit.

    """
    payload = {'api_key': api_key or DEFAULT_API_KEY}

    response = requests.get(NEOWS_URL + spk_id, params=payload)
    response.raise_for_status()

    orbital_data = response.json()['orbital_data']

    attractor = Sun
    a = float(orbital_data['semi_major_axis']) * u.AU
    ecc = float(orbital_data['eccentricity']) * u.one
    inc = float(orbital_data['inclination']) * u.deg
    raan = float(orbital_data['ascending_node_longitude']) * u.deg
    argp = float(orbital_data['perihelion_argument']) * u.deg
    m = float(orbital_data['mean_anomaly']) * u.deg
    nu = M_to_nu(m.to(u.rad), ecc)
    epoch = Time(float(orbital_data['epoch_osculation']), format='jd', scale='tdb')

    ss = Orbit.from_classical(attractor, a, ecc, inc,
                              raan, argp, nu, epoch, plane=Planes.EARTH_ECLIPTIC)
    return ss


def spk_id_from_name(name):
    """Return SPK-ID number given a small-body name.

    Retrieve and parse HTML from JPL Small Body Database
    to get SPK-ID.

    Parameters
    ----------
    name : str
        Small-body object name. Wildcards "*" and/or "?" can be used.

    Returns
    -------
    spk_id : str
        SPK-ID number.

    """
    payload = {'sstr': name, 'orb': '0', 'log': '0', 'old': '0', 'cov': '0', 'cad': '0'}

    response = requests.get(SBDB_URL, params=payload)
    response.raise_for_status()
    soup = BeautifulSoup(response.text, "html.parser")

    # page_identifier is used to check what type of response page we are working with.
    page_identifier = soup.find(attrs={"name": "top"})

    # If there is a 'table' sibling, the object was found.
    if page_identifier.find_next_sibling('table') is not None:
        data = page_identifier.find_next_sibling('table').table.find_all('td')

        complete_string = ''
        for string in data[1].stripped_strings:
            complete_string += string + ' '
        match = re.compile(r'Classification: ([\S\s]+) SPK-ID: (\d+)').match(complete_string)
        if match:
            return match.group(2)

    # If there is a 'center' sibling, it is a page with a list of possible objects
    elif page_identifier.find_next_sibling('center') is not None:
        object_list = page_identifier.find_next_sibling('center').table.find_all('td')
        bodies = ''
        obj_num = min(len(object_list), 3)
        for body in object_list[:obj_num]:
            bodies += body.string + '\n'
        raise ValueError(str(len(object_list)) + ' different bodies found:\n' + bodies)

    # If everything else failed
    raise ValueError('Object could not be found. You can visit: ' +
                     SBDB_URL + '?sstr=' + name + ' for more information.')


def orbit_from_name(name, api_key=None):
    """Return :py:class:`~poliastro.twobody.orbit.Orbit` given a name.

    Retrieve info from NASA NeoWS API, and therefore
    it only works with NEAs (Near Earth Asteroids).

    Parameters
    ----------
    name : str
        NEA name.
    api_key : str
        NASA OPEN APIs key (default: `DEMO_KEY`)

    Returns
    -------
    orbit : ~poliastro.twobody.orbit.Orbit
        NEA orbit.

    """
    spk_id = spk_id_from_name(name)
    if spk_id is not None:
        return orbit_from_spk_id(spk_id, api_key)
