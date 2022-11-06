---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.14.1
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Loading OMM and TLE satellite data

## How are satellite orbits disseminated?

![WTF is a TLE?](wtf-tle.jpg)

<em>(Image by <a href="https://twitter.com/fisadev/status/1309126214335037446/">@fisadev</a>)</em>

+++

- The only public source of orbital data for non-amateur Earth artificial satellites is the 18th Space Control Squadron (18 SPCS) from the United States Space Force (https://www.space-track.org/ and https://celestrak.org/NORAD/).
- The data is disseminated as _general perturbations_ (**GP**) orbital data following the _Simplified General Perturbations 4_ (**SGP4**) propagation model.
- Such GP orbital data has traditionally been transmitted in _Two-Line Element Sets_ or **TLEs**, although [in May 2020](https://celestrak.org/NORAD/documentation/gp-data-formats.php) it started being published in CCSDS _Orbit Mean-Elements Message_ (**OMM**) format as well.
- The advantage of SGP4 is that it's an _analytical model_ that includes the most important natural perturbations, and it's much faster than an equivalent _numerical model_. Therefore, it's good for initial orbit design based on low-precision requirements.

+++

However, it turns out that GP data in general, and TLEs in particular, are poorly understood even by professionals ([[1]](https://www.linkedin.com/posts/tom-johnson-32333a2_flawed-data-activity-6825845118990381056-yJX7), [[2]](https://twitter.com/flightclubio/status/1435303066085982209), [[3]](https://github.com/poliastro/poliastro/issues/1185)). The core issue is that TLEs and OMMs contain _Brouwer mean elements_, which **cannot be directly translated to osculating elements**.

From "Spacetrack Report #3":

> The **most important** point to be noted is that not just any prediction model will suffice. The NORAD element sets are "mean" values obtained by removing periodic variations in a particular way. In order to obtain good predictions, these periodic variations must be reconstructed (by the prediction model) in exactly the same way they were removed by NORAD.

From the Celestrak TLE FAQ:

> The elements in the two-line element sets are mean elements calculated to fit a set of observations using a specific model —the SGP4/SDP4 orbital model. Just as you shouldn't expect the arithmetic and geometric means of a set of data to have the same value, you shouldn't expect mean elements from different element sets —calculated using different orbital models— to have the same value. The short answer is that you cannot simply reformat the data unless you are willing to accept predictions with unpredictable errors.

From "Revisiting Spacetrack Report #3: Rev 3":

> Simply converting the orbital elements to an osculating state vector and propagating with a numerical propagator is equally invalid.

+++

Therefore, the **correct** way of using GP data is:

1. Propagate the TLE or OMM using the SGP4 algorithm, which produces cartesian elements (position and velocity) in the _True Equator Mean Equinox_ (TEME) reference frame.
2. Convert the resulting coordinates from TEME to the desired one (GCRS/ECI, ITRS/ECEF, other) (see note below).

<div class="alert alert-info">"Revisiting Spacetrack Report #3: Rev 3" warns that there is an inherent ambiguity in the reference frame used by the 18 SPCS to generate their GP data (TEME of Date vs TEME of Epoch) and there is no confirmation from official sources of which one is used. In any case, if you need sub-kilometer precision, GP data with SGP4 is not enough for you anyway.</div>

+++ {"slideshow": {"slide_type": "subslide"}}

# Loading general perturbations data

As explained in the [Orbit Mean-Elements Messages (OMMs) support assessment](https://opensatcom.org/2020/12/28/omm-assessment-sgp4-benchmarks/) deliverable of OpenSatCom, OMM input/output support in open source libraries is somewhat scattered. Luckily, [python-sgp4](https://pypi.org/project/sgp4/) supports reading OMM in CSV and XML format, as well as usual TLE and 3LE formats. On the other hand, Astropy has accurate transformations from TEME to other reference frames.

```{code-cell} ipython3
# From https://github.com/poliastro/poliastro/blob/main/contrib/satgpio.py
"""
Author: Juan Luis Cano Rodríguez

Code to read GP data from Celestrak using the HTTP API and python-sgp4.

Requires some extra dependencies:

  $ pip install httpx sgp4

This is similar to https://gitlab.com/librespacefoundation/python-satellitetle,
but uses the XML API instead and returns a `Satrec` object from sgp4 directly.

"""

import io
import json
import xml.etree.ElementTree as ET

import httpx
from sgp4 import exporter, omm
from sgp4.api import Satrec


def _generate_url(catalog_number, international_designator, name):
    params = {
        "CATNR": catalog_number,
        "INTDES": international_designator,
        "NAME": name,
    }
    param_names = [
        param_name
        for param_name, param_value in params.items()
        if param_value is not None
    ]
    if len(param_names) != 1:
        raise ValueError(
            "Specify exactly one of catalog_number, international_designator, or name"
        )
    param_name = param_names[0]
    param_value = params[param_name]
    url = (
        "https://celestrak.org/NORAD/elements/gp.php?"
        f"{param_name}={param_value}"
        "&FORMAT=XML"
    )
    return url


def _segments_from_query(url):
    response = httpx.get(url)
    response.raise_for_status()

    if response.text == "No GP data found":
        raise ValueError(
            f"Query '{url}' did not return any results, try a different one"
        )
    tree = ET.parse(io.StringIO(response.text))
    root = tree.getroot()

    yield from omm.parse_xml(io.StringIO(response.text))


def load_gp_from_celestrak(
    *, catalog_number=None, international_designator=None, name=None
):
    """Load general perturbations orbital data from Celestrak.

    Returns
    -------
    Satrec
        Orbital data from specified object.

    Notes
    -----
    This uses the OMM XML format from Celestrak as described in [1]_.

    References
    ----------
    .. [1] Kelso, T.S. "A New Way to Obtain GP Data (aka TLEs)"
       https://celestrak.org/NORAD/documentation/gp-data-formats.php

    """
    # Assemble query, raise an error if malformed
    url = _generate_url(catalog_number, international_designator, name)

    # Make API call, raise an error if data is malformed
    for segment in _segments_from_query(url):
        # Initialize and return Satrec object
        sat = Satrec()
        omm.initialize(sat, segment)

        yield sat


def print_sat(sat, name):
    """Prints Satrec object in convenient form."""
    print(json.dumps(exporter.export_omm(sat, name), indent=2))
```

```{code-cell} ipython3
sat = list(load_gp_from_celestrak(name="ISS (Zarya)"))[0]
print_sat(sat, "ISS (Zarya)")
```

The generator `load_gp_from_celestrak` might yield more than one satellite:

```{code-cell} ipython3
len(list(load_gp_from_celestrak(name="COSMOS 1408 DEB")))
```

# Creating ephemerides from general perturbations data

Now that we know how to load GP data, let's propagate it using the SGP4 algorithm:

```{code-cell} ipython3
from astropy import units as u
from astropy.time import Time
```

```{code-cell} ipython3
now = Time.now()
now.jd1, now.jd2
```

```{code-cell} ipython3
error, r, v = sat.sgp4(now.jd1, now.jd2)
assert error == 0
```

```{code-cell} ipython3
r << u.km
```

```{code-cell} ipython3
v << (u.km / u.s)
```

Note that `Satrec` also supports propagating for arrays of time values using the `.sgp4_array` method:

```{code-cell} ipython3
import numpy as np

from astropy.coordinates import CartesianRepresentation, CartesianDifferential

from poliastro.util import time_range
```

```{code-cell} ipython3
times = time_range(now, end=now + (1 << u.h), num_values=3)
```

```{code-cell} ipython3
errors, rs, vs = sat.sgp4_array(times.jd1, times.jd2)
assert (errors == 0).all()
```

```{code-cell} ipython3
CartesianRepresentation(rs << u.km, xyz_axis=-1)
```

```{code-cell} ipython3
CartesianDifferential(vs << (u.km / u.s), xyz_axis=-1)
```

Mixing all together, and leveraging poliastro `Ephem` objects, we can compute ephemerides from GP orbital data:

```{code-cell} ipython3
from warnings import warn

from astropy.coordinates import TEME, GCRS

from poliastro.ephem import Ephem
from poliastro.frames import Planes


def ephem_from_gp(sat, times):
    errors, rs, vs = sat.sgp4_array(times.jd1, times.jd2)
    if not (errors == 0).all():
        warn(
            "Some objects could not be propagated, "
            "proceeding with the rest",
            stacklevel=2,
        )
        rs = rs[errors == 0]
        vs = vs[errors == 0]
        times = times[errors == 0]

    cart_teme = CartesianRepresentation(
        rs << u.km,
        xyz_axis=-1,
        differentials=CartesianDifferential(
            vs << (u.km / u.s),
            xyz_axis=-1,
        ),
    )
    cart_gcrs = (
        TEME(cart_teme, obstime=times)
        .transform_to(GCRS(obstime=times))
        .cartesian
    )

    return Ephem(cart_gcrs, times, plane=Planes.EARTH_EQUATOR)
```

```{code-cell} ipython3
iss_ephem = ephem_from_gp(sat, time_range(now, end=now + (3 << u.h)))
iss_ephem
```

And plot it!

```{code-cell} ipython3
from poliastro.bodies import Earth
from poliastro.plotting import OrbitPlotter

plotter = OrbitPlotter(backend_name="plotly3D")
plotter.set_attractor(Earth)
plotter.plot_ephem(iss_ephem, color="#333", label="ISS", trail=True)

plotter.show()
```

We can also plot the trajectory of the ISS for 90 minutes, plus the trajectories of the first 25 debris fragments:

```{code-cell} ipython3
from itertools import islice

epochs = time_range(now, end=now + (90 << u.minute))

iss_ephem = ephem_from_gp(sat, epochs)

plotter = OrbitPlotter(backend_name="plotly3D")
plotter.set_attractor(Earth)
plotter.plot_ephem(iss_ephem, color="#333", label="ISS", trail=True)

for debris_fragment in islice(
    load_gp_from_celestrak(name="COSMOS 1408 DEB"), 25
):
    debris_ephem = ephem_from_gp(debris_fragment, epochs)
    plotter.plot_ephem(
        debris_ephem, color="#666", label=debris_fragment.satnum, trail=True
    )

plotter.show()
```
