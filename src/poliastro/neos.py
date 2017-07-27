'''Code related to NEOs

All new functions coded as part of SOCIS 2017 proposal

'''
import re
from bs4 import BeautifulSoup
import requests

import astropy.units as u

from poliastro.twobody.orbit import Orbit
from poliastro.bodies import Sun
from poliastro.twobody.angles import M_to_nu

def get_orbit_from_spk_id(spk_id, epoch=None):
    """Return `~poliastro.twobody.orbit.Orbit` given a SPK-ID.

    Retrieve info from NASA NeoWS API, and therefore
    it only works with NEAs (Near Earth Asteroids)

    Parameters
    ----------
    spk_id : str
        SPK-ID number, which is given to each body by JPL.

    Returns
    -------
    orbit : ~poliastro.twobody.orbit.Orbit
        NEA orbit..

    """
    url = 'https://api.nasa.gov/neo/rest/v1/neo/'
    payload = {'api_key' : 'DEMO_KEY'}
    try:
        response = requests.get(url + spk_id, params=payload)
        if response.status_code == 200:
            orbital_data = response.json()['orbital_data']

            attractor = Sun
            a = float(orbital_data['semi_major_axis']) * u.AU
            ecc = float(orbital_data['eccentricity']) * u.one
            inc = float(orbital_data['inclination']) * u.deg
            anl = float(orbital_data['ascending_node_longitude']) * u.deg
            parg = float(orbital_data['perihelion_argument']) * u.deg
            M = float(orbital_data['mean_anomaly']) * u.deg
            nu = M_to_nu(M.to(u.rad), ecc)

            if epoch is None:
                return Orbit.from_classical(attractor, a, ecc, inc,
                                        anl, parg, nu)
            else:
                return Orbit.from_classical(attractor, a, ecc, inc,
                                        anl, parg, nu, epoch)
        elif response.status_code == 404:
            print('\nError 404. Object ID does not exists, or it does not match a NEA.')
        else:
            print('\nForbidden or unauthorized request: ', response.status_code)
        return None
    except requests.exceptions.RequestException:
        print('\nSome kind of connection error ocurred, please try later.')
        return None
    except KeyError:
        print('Internal error, please report it to Poliastro team')
        return None

def get_spk_id_from_name(name):
    '''Return SPK-ID number given a small-body name.
    
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

    '''
    url = 'https://ssd.jpl.nasa.gov/sbdb.cgi'
    headers = {'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_1) \
        AppleWebKit/537.36 (KHTML, like Gecko) Chrome/39.0.2171.95 Safari/537.36'}
    payload = {'sstr' : name, 'orb' : '0', 'log' : '0', 'old' : '0', 'cov' : '0', 'cad' : '0'}
    try:
        response = requests.get(url, params=payload, headers=headers)
        soup = BeautifulSoup(response.text, "html.parser")

        # page_identifier is used to check what type of response page we are working with.
        page_identifier = soup.find(attrs={"name": "top"})

        # If there is a 'table' sibling, the object was found.
        if page_identifier.find_next_sibling('table') is not None:
            data = page_identifier.find_next_sibling('table').table.find_all('td')

            complete_string = ''
            for string in data[1].stripped_strings:
                complete_string += string + ' '
            regex = re.compile('Classification: ([\S\s]+) SPK-ID: (\d+)')
            return regex.match(complete_string).group(2)
        # If there is a 'center' sibling, it is a page with a list of possible objects
        elif page_identifier.find_next_sibling('center') is not None:
            object_list = page_identifier.find_next_sibling('center').table.find_all('td')
            print('The following ', len(object_list), ' objects match your search string:\n')
            for body in object_list:
                #body.a['href']
                print(body.string)
            print('\nPlease select one of them')
        else:
            print('Object could not be found. You can visit: ',
                 url + "?sstr=" + name, " for more information.")
        return None
    except requests.exceptions.RequestException:
        print('Some kind of connection error ocurred, please try later.')
        return None
    except (AttributeError, KeyError):
        print('Object could not be found. You can visit: ',
             url + "?sstr=" + name, " for more information.")
        return None

def get_orbit_from_name(name):
    '''Return `~poliastro.twobody.orbit.Orbit` given a name.

    Retrieve info from NASA NeoWS API, and therefore
    it only works with NEAs (Near Earth Asteroids)

    Parameters
    ----------
    name : str
        NEA name.

    Returns
    -------
    orbit : ~poliastro.twobody.orbit.Orbit
        NEA orbit.

    '''
    return get_orbit_from_spk_id(get_spk_id_from_name(name))
