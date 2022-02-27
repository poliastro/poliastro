#!/usr/bin/env python
# -*- coding: utf-8 -*-

from gzip import decompress
import json
import os
from typing import Union
from urllib.request import urlopen, Request

from astropy.time import Time
import astropy.units as u
import numpy as np

from poliastro.bodies import Sun
from poliastro.frames import Planes
from poliastro.twobody import Orbit
from poliastro.core.angles import (
    D_to_nu,
    E_to_nu,
    F_to_nu,
    M_to_D,
    M_to_E,
    M_to_F,
)

CWD = os.path.dirname(__file__)
NEA_MPC_URL = 'https://minorplanetcenter.net/Extended_Files/neam00_extended.json.gz'
DATA_FN = os.path.join(CWD, 'data.json')

def benchmark():

    orbs = [_orbit_from_mpc(body) for body in _get_data()[:10]] # TODO limit

    print(len(orbs), type(orbs))

def _download(down_url: str, binary: bool = True) -> Union[str, bytes]:

    httprequest = Request(down_url)

    with urlopen(httprequest) as response:
        assert response.status == 200
        data = response.read()

    if not binary:
        return data.decode('utf-8')

    return data

def _get_data() -> list[dict]:

    if not os.path.exists(DATA_FN):
        raw = _download(NEA_MPC_URL)
        with open(DATA_FN, mode = 'wb') as f:
            f.write(raw)
    else:
        with open(DATA_FN, mode = 'rb') as f:
            raw = f.read()

    return json.loads(decompress(raw))

def _orbit_from_mpc(body: dict) -> Orbit:

    nu = _true_anomaly_from_mean(
        ecc = body['e'],
        M = float((body['M'] * u.deg).to(u.rad).value),
    ) * u.rad

    if not -np.pi * u.rad <= nu < np.pi * u.rad:
        nu = ((nu + np.pi * u.rad) % (2 * np.pi * u.rad) - np.pi * u.rad).to(
            nu.unit
        )

    return Orbit.from_classical(
        Sun,
        a = (body['a'] * u.AU).to(u.m),
        ecc = body['e'] * u.one,
        inc = (body['i'] * u.deg).to(u.rad),
        raan = (body['Node'] * u.deg).to(u.rad),
        argp = (body['Peri'] * u.deg).to(u.rad),
        nu = nu,
        epoch = Time(body["Epoch"], format = 'jd'),
        plane = Planes.EARTH_ECLIPTIC,
    )

def _true_anomaly_from_mean(ecc: float, M: float) -> float:

    if ecc < 1:
        M = (M + np.pi) % (2 * np.pi) - np.pi
        return E_to_nu(M_to_E(M, ecc), ecc)
    elif ecc == 1:
        return D_to_nu(M_to_D(M))
    else:
        return F_to_nu(M_to_F(M, ecc), ecc)

if __name__ == '__main__':

    benchmark()
