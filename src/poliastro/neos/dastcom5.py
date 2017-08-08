"""NEOs orbit from DASTCOM5 database.

"""
import os
import re
import zipfile
import urllib.request

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.time import Time

from poliastro.twobody.orbit import Orbit
from poliastro.bodies import Sun
from poliastro.twobody.angles import M_to_nu

AST_DTYPE = np.dtype([
    ("NO", np.int32),
    ("NOBS", np.int32),
    ("OBSFRST", np.int32),
    ("OBSLAST", np.int32),
    ("EPOCH", np.float64),
    ("CALEPO", np.float64),
    ("MA", np.float64),
    ("W", np.float64),
    ("OM", np.float64),
    ("IN", np.float64),
    ("EC", np.float64),
    ("A", np.float64),
    ("QR", np.float64),
    ("TP", np.float64),
    ("TPCAL", np.float64),
    ("TPFRAC", np.float64),
    ("SOLDAT", np.float64),
    ("SRC1", np.float64),
    ("SRC2", np.float64),
    ("SRC3", np.float64),
    ("SRC4", np.float64),
    ("SRC5", np.float64),
    ("SRC6", np.float64),
    ("SRC7", np.float64),
    ("SRC8", np.float64),
    ("SRC9", np.float64),
    ("SRC10", np.float64),
    ("SRC11", np.float64),
    ("SRC12", np.float64),
    ("SRC13", np.float64),
    ("SRC14", np.float64),
    ("SRC15", np.float64),
    ("SRC16", np.float64),
    ("SRC17", np.float64),
    ("SRC18", np.float64),
    ("SRC19", np.float64),
    ("SRC20", np.float64),
    ("SRC21", np.float64),
    ("SRC22", np.float64),
    ("SRC23", np.float64),
    ("SRC24", np.float64),
    ("SRC25", np.float64),
    ("SRC26", np.float64),
    ("SRC27", np.float64),
    ("SRC28", np.float64),
    ("SRC29", np.float64),
    ("SRC30", np.float64),
    ("SRC31", np.float64),
    ("SRC32", np.float64),
    ("SRC33", np.float64),
    ("SRC34", np.float64),
    ("SRC35", np.float64),
    ("SRC36", np.float64),
    ("SRC37", np.float64),
    ("SRC38", np.float64),
    ("SRC39", np.float64),
    ("SRC40", np.float64),
    ("SRC41", np.float64),
    ("SRC42", np.float64),
    ("SRC43", np.float64),
    ("SRC44", np.float64),
    ("SRC45", np.float64),
    ("PRELTV", np.int8),
    ("SPHMX3", np.int8),
    ("SPHMX5", np.int8),
    ("JGSEP", np.int8),
    ("TWOBOD", np.int8),
    ("NSATS", np.int8),
    ("UPARM", np.int8),
    ("LSRC", np.int8),
    ("NDEL", np.int16),
    ("NDOP", np.int16),
    ("H", np.float32),
    ("G", np.float32),
    ("A1", np.float32),
    ("A2", np.float32),
    ("A3", np.float32),
    ("R0", np.float32),
    ("ALN", np.float32),
    ("NM", np.float32),
    ("NN", np.float32),
    ("NK", np.float32),
    ("LGK", np.float32),
    ("RHO", np.float32),
    ("AMRAT", np.float32),
    ("ALF", np.float32),
    ("DEL", np.float32),
    ("SPHLM3", np.float32),
    ("SPHLM5", np.float32),
    ("RP", np.float32),
    ("GM", np.float32),
    ("RAD", np.float32),
    ("EXTNT1", np.float32),
    ("EXTNT2", np.float32),
    ("EXTNT3", np.float32),
    ("MOID", np.float32),
    ("ALBEDO", np.float32),
    ("BVCI", np.float32),
    ("UBCI", np.float32),
    ("IRCI", np.float32),
    ("RMSW", np.float32),
    ("RMSU", np.float32),
    ("RMSN", np.float32),
    ("RMSNT", np.float32),
    ("RMSH", np.float32),
    ("EQUNOX", "|S4"),
    ("PENAM", "|S6"),
    ("SBNAM", "|S12"),
    ("SPTYPT", "|S5"),
    ("SPTYPS", "|S5"),
    ("DARC", "|S9"),
    ("COMNT1", "|S41"),
    ("COMNT2", "|S80"),
    ("DESIG", "|S13"),
    ("ASTEST", "|S8"),
    ("IREF", "|S10"),
    ("ASTNAM", "|S18"),
])

COM_DTYPE = np.dtype([
    ("NO", np.int32),
    ("NOBS", np.int32),
    ("OBSFRST", np.int32),
    ("OBSLAST", np.int32),
    ("EPOCH", np.float64),
    ("CALEPO", np.float64),
    ("MA", np.float64),
    ("W", np.float64),
    ("OM", np.float64),
    ("IN", np.float64),
    ("EC", np.float64),
    ("A", np.float64),
    ("QR", np.float64),
    ("TP", np.float64),
    ("TPCAL", np.float64),
    ("TPFRAC", np.float64),
    ("SOLDAT", np.float64),
    ("SRC1", np.float64),
    ("SRC2", np.float64),
    ("SRC3", np.float64),
    ("SRC4", np.float64),
    ("SRC5", np.float64),
    ("SRC6", np.float64),
    ("SRC7", np.float64),
    ("SRC8", np.float64),
    ("SRC9", np.float64),
    ("SRC10", np.float64),
    ("SRC11", np.float64),
    ("SRC12", np.float64),
    ("SRC13", np.float64),
    ("SRC14", np.float64),
    ("SRC15", np.float64),
    ("SRC16", np.float64),
    ("SRC17", np.float64),
    ("SRC18", np.float64),
    ("SRC19", np.float64),
    ("SRC20", np.float64),
    ("SRC21", np.float64),
    ("SRC22", np.float64),
    ("SRC23", np.float64),
    ("SRC24", np.float64),
    ("SRC25", np.float64),
    ("SRC26", np.float64),
    ("SRC27", np.float64),
    ("SRC28", np.float64),
    ("SRC29", np.float64),
    ("SRC30", np.float64),
    ("SRC31", np.float64),
    ("SRC32", np.float64),
    ("SRC33", np.float64),
    ("SRC34", np.float64),
    ("SRC35", np.float64),
    ("SRC36", np.float64),
    ("SRC37", np.float64),
    ("SRC38", np.float64),
    ("SRC39", np.float64),
    ("SRC40", np.float64),
    ("SRC41", np.float64),
    ("SRC42", np.float64),
    ("SRC43", np.float64),
    ("SRC44", np.float64),
    ("SRC45", np.float64),
    ("SRC46", np.float64),
    ("SRC47", np.float64),
    ("SRC48", np.float64),
    ("SRC49", np.float64),
    ("SRC50", np.float64),
    ("SRC51", np.float64),
    ("SRC52", np.float64),
    ("SRC53", np.float64),
    ("SRC54", np.float64),
    ("SRC55", np.float64),
    ("PRELTV", np.int8),
    ("SPHMX3", np.int8),
    ("SPHMX5", np.int8),
    ("JGSEP", np.int8),
    ("TWOBOD", np.int8),
    ("NSATS", np.int8),
    ("UPARM", np.int8),
    ("LSRC", np.int8),
    ("IPYR", np.int16),
    ("NDEL", np.int16),
    ("NDOP", np.int16),
    ("NOBSMT", np.int16),
    ("NOBSMN", np.int16),
    ("H", np.float32),
    ("G", np.float32),
    ("M1 (MT)", np.float32),
    ("M2 (MN)", np.float32),
    ("K1 (MTSMT)", np.float32),
    ("K2 (MNSMT)", np.float32),
    ("PHCOF (MNP)", np.float32),
    ("A1", np.float32),
    ("A2", np.float32),
    ("A3", np.float32),
    ("DT", np.float32),
    ("R0", np.float32),
    ("ALN", np.float32),
    ("NM", np.float32),
    ("NN", np.float32),
    ("NK", np.float32),
    ("S0", np.float32),
    ("TCL", np.float32),
    ("RHO", np.float32),
    ("AMRAT", np.float32),
    ("AJ1", np.float32),
    ("AJ2", np.float32),
    ("ET1", np.float32),
    ("ET2", np.float32),
    ("DTH", np.float32),
    ("ALF", np.float32),
    ("DEL", np.float32),
    ("SPHLM3", np.float32),
    ("SPHLM5", np.float32),
    ("RP", np.float32),
    ("GM", np.float32),
    ("RAD", np.float32),
    ("EXTNT1", np.float32),
    ("EXTNT2", np.float32),
    ("EXTNT3", np.float32),
    ("MOID", np.float32),
    ("ALBEDO", np.float32),
    ("RMSW", np.float32),
    ("RMSU", np.float32),
    ("RMSN", np.float32),
    ("RMSNT", np.float32),
    ("RMSMT", np.float32),
    ("RMSMN", np.float32),
    ("EQUNOX", "|S4"),
    ("PENAM", "|S6"),
    ("SBNAM", "|S12"),
    ("DARC", "|S9"),
    ("COMNT3", "|S49"),
    ("COMNT2", "|S80"),
    ("DESIG", "|S13"),
    ("COMEST", "|S14"),
    ("IREF", "|S10"),
    ("COMNAM", "|S29"),
])

POLIASTRO_LOCAL_PATH = os.path.join(os.path.expanduser('~'), '.poliastro')
DBS_LOCAL_PATH = os.path.join(POLIASTRO_LOCAL_PATH, 'dastcom5', 'dat')
AST_DB_PATH = os.path.join(DBS_LOCAL_PATH, 'dast5_le.dat')
COM_DB_PATH = os.path.join(DBS_LOCAL_PATH, 'dcom5_le.dat')

FTP_DB_URL = 'ftp://ssd.jpl.nasa.gov/pub/ssd/'


def asteroid_db_from_path(filepath=AST_DB_PATH):
    """Return complete DASTCOM5 asteroid database.

    Parameters
    ----------
    filepath : path-like object
        dast5_le.dat file path.

    Returns
    -------
    database : numpy.ndarray
        Database with custom dtype.

    """
    with open(filepath, "rb") as f:
        f.seek(835, os.SEEK_SET)
        data = np.fromfile(f, dtype=AST_DTYPE)
    return data


def comet_db_from_path(filepath=COM_DB_PATH):
    """Return complete DASTCOM5 comet database.

    Parameters
    ----------
    filepath : path-like object
        dcom5_le.dat file path.

    Returns
    -------
    database : numpy.ndarray
        Database with custom dtype.

    """
    with open(filepath, "rb") as f:
        f.seek(976, os.SEEK_SET)
        data = np.fromfile(f, dtype=COM_DTYPE)
    return data


def orbit_from_name(name):
    """Return :py:class:`~poliastro.twobody.orbit.Orbit` given a name.

    Retrieve info from JPL DASTCOM5 database

    Parameters
    ----------
    name : str
        NEO name.

    Returns
    -------
    orbit : ~poliastro.twobody.orbit.Orbit
        NEO orbit.

    """
    records = record_from_name(name)
    orbits = []

    for record in records:
        body_data = read_record(record)
        a = body_data['A'].item() * u.au
        ecc = body_data['EC'].item() * u.one
        inc = body_data['IN'].item() * u.deg
        raan = body_data['OM'].item() * u.deg
        argp = body_data['W'].item() * u.deg
        m = body_data['MA'].item() * u.deg
        nu = M_to_nu(m, ecc)
        epoch = Time(body_data['EPOCH'].item(), format='jd')

        orbits.append(Orbit.from_classical(Sun, a, ecc, inc, raan, argp, nu, epoch))
    return orbits


def record_from_name(name, dirpath=DBS_LOCAL_PATH):
    """Search `dastcom.idx` and return logical records that match a given string.

    (Body name, SPK-ID, or alternative designations
         can be used).

    Parameters
    ----------
    name : str
        Body name.
    dirpath : path-like object
        dastcom.idx dir.

    Returns
    -------
    records : list (int)
        DASTCOM5 database logical records matching str.

    """
    idx_path = os.path.join(dirpath, 'dastcom.idx')
    records = []
    with open(idx_path, 'r') as inF:
        for line in inF:
            if re.search(r"\b" + name.casefold() + r"\b", line.casefold()):
                records.append(int(line[:6].lstrip()))
    return records


def read_headers(dirpath=DBS_LOCAL_PATH):
    """Read `DASTCOM5` headers and return asteroid and comet headers.

    Headers are two numpy arrays with custom dtype.

    Parameters
    ----------
    dirpath : path-like object
        DASTCOM5 .dat files dir.

    Returns
    -------
    ast_header, com_header : tuple (numpy.ndarray)
        DASTCOM5 headers.

    """

    ast_path = os.path.join(dirpath, 'dast5_le.dat')
    ast_dtype = np.dtype([
        ("IBIAS1", np.int32),
        ("BEGINP1", "|S8"),
        ("BEGINP2", "|S8"),
        ("BEGINP3", "|S8"),
        ("ENDPT1", "|S8"),
        ("ENDPT2", "|S8"),
        ("ENDPT3", "|S8"),
        ("CALDATE", "|S19"),
        ("JDDATE", np.float64),
        ("FTYP", "|S1"),
        ("BYTE2A", np.int16),
        ("IBIAS0", np.int32),
    ])

    with open(ast_path, 'rb') as f:
        ast_header = np.fromfile(f, dtype=ast_dtype, count=1)

    com_path = os.path.join(dirpath, 'dcom5_le.dat')
    com_dtype = np.dtype([
        ("IBIAS2", np.int32),
        ("BEGINP1", "|S8"),
        ("BEGINP2", "|S8"),
        ("BEGINP3", "|S8"),
        ("ENDPT1", "|S8"),
        ("ENDPT2", "|S8"),
        ("ENDPT3", "|S8"),
        ("CALDATE", "|S19"),
        ("JDDATE", np.float64),
        ("FTYP", "|S1"),
        ("BYTE2C", np.int16),
    ])

    with open(com_path, 'rb') as f:
        com_header = np.fromfile(f, dtype=com_dtype, count=1)

    return ast_header, com_header


def read_record(record, dirpath=DBS_LOCAL_PATH):
    """Read `DASTCOM5` record and return body data.

    Body data consists of numpy array with custom dtype.

    Parameters
    ----------
    record : int
        Body name.
    dirpath : path-like object
        dastcom.idx dir.

    Returns
    -------
    body_data : numpy.ndarray
        Body information.

    """
    ast_header, com_header = read_headers()
    ast_path = os.path.join(dirpath, 'dast5_le.dat')
    com_path = os.path.join(dirpath, 'dcom5_le.dat')
    # ENDPT1 indicates end of numbered asteroids records
    if record <= int(ast_header['ENDPT2'][0].item()):
        # ENDPT2 indicates end of unnumbered asteroids records
        if record <= int(ast_header['ENDPT1'][0].item()):
            # phis_rec = record_size * (record_number - IBIAS - 1 (header record))
            phis_rec = 835 * (record - ast_header['IBIAS0'][0].item() - 1)
        else:
            phis_rec = 835 * (record - ast_header['IBIAS1'][0].item() - 1)

        with open(ast_path, "rb") as f:
            f.seek(phis_rec, os.SEEK_SET)
            body_data = np.fromfile(f, dtype=AST_DTYPE, count=1)
    else:
        phis_rec = 976 * (record - com_header['IBIAS2'][0].item() - 1)
        with open(com_path, "rb") as f:
            f.seek(phis_rec, os.SEEK_SET)
            body_data = np.fromfile(f, dtype=COM_DTYPE, count=1)
    return body_data


def download_dastcom5(path=POLIASTRO_LOCAL_PATH):
    """Downloads DASTCOM5 database to specified path.

    The function downloads and unzip DASTCOM5 file in a given path.

    Parameters
    ----------
    path : path-like object
        DASTCOM5 dir.

    """

    dastcom5_dir = os.path.join(path, 'dastcom5')
    dastcom5_zip_path = os.path.join(path, 'dastcom5.zip')

    if os.path.isdir(dastcom5_dir):
        raise FileExistsError('dastcom5 is already created in ' + os.path.abspath(dastcom5_dir))
    if not zipfile.is_zipfile(dastcom5_zip_path):
        if not os.path.isdir(path):
            os.makedirs(path)

        urllib.request.urlretrieve(FTP_DB_URL + 'dastcom5.zip',
                                   dastcom5_zip_path, _show_download_progress)
    with zipfile.ZipFile(dastcom5_zip_path) as myzip:
        myzip.extractall(path)


def _show_download_progress(transferred, block, totalsize):
    trans_mb = transferred*block/(1024*1024)
    total_mb = totalsize/(1024*1024)
    print('%.2f MB / %.2f MB' % (trans_mb, total_mb), end='\r', flush=True)


def entire_db_from_path(path=DBS_LOCAL_PATH):
    """Return complete DASTCOM5 database.

    Functions merge asteroid and comet database, only with fields
    related to orbital data, discarding the rest.

    Parameters
    ----------
    path : path-like object
        DASTCOM5 /dat folder path.

    Returns
    -------
    database : numpy.ndarray
        Database with custom dtype.

    """
    ast_path = os.path.join(path, 'dast5_le.dat')
    com_path = os.path.join(path, 'dcom5_le.dat')

    ast_database = asteroid_db_from_path(ast_path)
    com_database = comet_db_from_path(com_path)

    ast_database = pd.DataFrame(ast_database[list(ast_database.dtype.names[:17])
                                             + list(ast_database.dtype.names[-4:-3])
                                             + list(ast_database.dtype.names[-2:])])
    ast_database.rename(columns={'ASTNAM': 'NAME', 'NO': 'NUMBER', 'CALEPO': 'CALEPOCH'}, inplace=True)
    com_database = pd.DataFrame(com_database[list(com_database.dtype.names[:17])
                                             + list(com_database.dtype.names[-4:-3])
                                             + list(com_database.dtype.names[-2:])])
    com_database.rename(columns={'COMNAM': 'NAME', 'NO': 'NUMBER', 'CALEPO': 'CALEPOCH'}, inplace=True)
    df = ast_database.append(com_database, ignore_index=True)
    df[['NAME', 'DESIG', 'IREF']] = df[['NAME', 'DESIG', 'IREF']].apply(lambda x: x.str.strip().str.decode('utf-8'))
    return df
