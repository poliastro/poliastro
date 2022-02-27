#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gc
from gzip import decompress
import json
from multiprocessing import Pool, cpu_count
import os
from time import time
from typing import Union
from urllib.request import urlopen, Request

from astropy.time import Time
import astropy.units as u
import numpy as np
import plotly.graph_objects as go

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
MPC_URL = 'https://minorplanetcenter.net/Extended_Files/neam00_extended.json.gz'
DATA_FN = os.path.join(CWD, 'data.json')
PLOT_FN = os.path.join(CWD, 'benchmark.htm')


def benchmark(
    limit_exp: float = 4.0, # 10 ** limit_exp
    limit_exp_step: float = 0.5,
    plot_fn: str = PLOT_FN,
    data_fn: str = DATA_FN,
    data_url: str = MPC_URL,
):
    "run and plot simple benchmark"

    assert limit_exp >= 0
    steps = (10 ** np.arange(0, limit_exp + limit_exp_step, limit_exp_step)).astype('i8')

    bodies = _get_bodies(data_fn = data_fn, data_url = data_url)
    limit = int(steps[-1])
    assert limit <= len(bodies)
    orbs = _orbs_from_mpc(bodies[:limit])

    epoch = Time('2022-01-01 00:00')

    _ = [orb.propagate(epoch) for orb in orbs] # JIT warm-up on all code-paths

    gc_flag = gc.isenabled()
    if gc_flag:
        gc.disable() # consistent benchmarks

    times = []
    for idx, iterations in enumerate(steps):

        print(f'Benchmark {idx:d} of {len(steps)}, {iterations:d} iteratations ...')
        selected_orbs = orbs[:int(iterations)]

        start_time = time() # START

        _ = [orb.propagate(epoch) for orb in selected_orbs]

        times.append(time() - start_time) # STOP

        gc.collect()
        print(f'... finished in {times[-1]:e} seconds ({times[-1]/iterations:e} seconds per iteration).')

    if gc_flag:
        gc.enable()

    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x = steps,
            y = np.array(times) / steps,
            name = 'Iterative',
            line = dict(width = 2)
        )
    )

    fig.update_layout(
        title = 'Simple Orbit Propagation Benchmark',
        xaxis_title = 'Number of Propagations',
        xaxis_type = "log",
        yaxis_title = 'Time per Propagation [s]',
        yaxis_type = "log",
    )

    fig.write_html(plot_fn, auto_open = True)


def _download(url: str, binary: bool = True) -> Union[str, bytes]:
    "simple, dependency-free http fetch"

    httprequest = Request(url)

    with urlopen(httprequest) as response:
        assert response.status == 200
        data = response.read()

    if binary:
        return data

    return data.decode('utf-8')


def _get_bodies(data_fn: str, data_url: str) -> list[dict]:
    "load MPC data from disk and/or internet"

    if not os.path.exists(data_fn):
        raw = _download(data_url)
        with open(data_fn, mode = 'wb') as f:
            f.write(raw)
    else:
        with open(data_fn, mode = 'rb') as f:
            raw = f.read()

    return json.loads(decompress(raw))


def _orb_from_mpc(body: dict) -> Orbit:
    "convert MPC single body dict into Orbit object"

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


def _orbs_from_mpc(bodies: list[dict]) -> list[Orbit]:
    "convert list of MPC single body dicts into list of Orbit objects (parallel)"

    with Pool(cpu_count()) as pool:
        tasks = [
            pool.apply_async(
                func = _orb_from_mpc,
                args = (body,),
            )
            for body in bodies
        ]
        orbs = [task.get() for task in tasks]

    return orbs


def _true_anomaly_from_mean(ecc: float, M: float) -> float:
    "see issue #1013"

    if ecc < 1:
        M = (M + np.pi) % (2 * np.pi) - np.pi
        return E_to_nu(M_to_E(M, ecc), ecc)
    if ecc == 1:
        return D_to_nu(M_to_D(M))
    return F_to_nu(M_to_F(M, ecc), ecc)


if __name__ == '__main__':

    benchmark()
