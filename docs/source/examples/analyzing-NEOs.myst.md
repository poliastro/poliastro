---
jupytext:
  encoding: '# -*- coding: utf-8 -*-'
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.14.4
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Analyzing NEOs

+++

NEO stands for near-Earth object. The Center for NEO Studies ([CNEOS](http://cneos.jpl.nasa.gov/)) defines NEOs as comets and asteroids that have been nudged by the gravitational attraction of nearby planets into orbits that allow them to enter the Earth’s neighborhood.

And what does "near" exactly mean? In terms of orbital elements, asteroids and comets can be considered NEOs if their perihelion (orbit point which is nearest to the Sun) is less than 1.3 au = 1.945 * 10<sup>8</sup> km from the Sun.

```{code-cell} ipython3
from astropy import time

from poliastro.bodies import Earth
from poliastro.frames import Planes
from poliastro.plotting import OrbitPlotter
from poliastro.twobody.orbit import Orbit
```

## Small Body Database (SBDB)

```{code-cell} ipython3
eros = Orbit.from_sbdb("Eros")
eros.plot(label="Eros")
```

You can also search by IAU number or SPK-ID (there is a faster `neows.orbit_from_spk_id()` function in that case, although):

```{code-cell} ipython3
:tags: [nbsphinx-thumbnail]

ganymed = Orbit.from_sbdb("1036")  # Ganymed IAU number
amor = Orbit.from_sbdb("2001221")  # Amor SPK-ID
eros = Orbit.from_sbdb("2000433")  # Eros SPK-ID

frame = OrbitPlotter(plane=Planes.EARTH_ECLIPTIC)
frame.plot(ganymed, label="Ganymed")
frame.plot(amor, label="Amor")
frame.plot(eros, label="Eros")
```

You can use the wildcards from that browser: `*` and `?`.

+++

<div class="alert alert-info">Keep it in mind that `from_sbdb()` can only return one Orbit, so if several objects are found with that name, it will raise an error with the different bodies:</div>

```{code-cell} ipython3
try:
    Orbit.from_sbdb("*alley")
except ValueError as err:
    print(err)
```

<div class="alert alert-info">Note that the epoch is provided by the service itself, so if you need orbit on another epoch, you have to propagate it:</div>

```{code-cell} ipython3
eros.epoch.iso
```

```{code-cell} ipython3
epoch = time.Time(2458000.0, scale="tdb", format="jd")
eros_november = eros.propagate(epoch)
eros_november.epoch.iso
```
