#!/usr/bin/env python
# -*- coding: utf-8 -*-

from gzip import decompress
import json
import os
from typing import Union
from urllib.request import urlopen, Request

CWD = os.path.dirname(__file__)
NEA_MPC_URL = 'https://minorplanetcenter.net/Extended_Files/neam00_extended.json.gz'
DATA_FN = os.path.join(CWD, 'data.json')

def benchmark():

    data = _get_data()

    print(len(data), type(data))

def _download(down_url: str, mode: str = "binary") -> Union[str, bytes]:

    assert mode in ("text", "binary")
    assert isinstance(down_url, str)

    httprequest = Request(down_url)

    with urlopen(httprequest) as response:
        assert response.status == 200
        data = response.read()

    if mode == 'text':
        return data.decode('utf-8')

    return data # mode == 'binary'

def _get_data() -> list[dict]:

    if not os.path.exists(DATA_FN):
        raw = _download(NEA_MPC_URL)
        with open(DATA_FN, mode = 'wb') as f:
            f.write(raw)
    else:
        with open(DATA_FN, mode = 'rb') as f:
            raw = f.read()

    return json.loads(decompress(raw))

if __name__ == '__main__':

    benchmark()
