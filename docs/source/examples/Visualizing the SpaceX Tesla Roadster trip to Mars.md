---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.5.0
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Visualizing the SpaceX Tesla Roadster trip to Mars

```python
from astropy.time import Time
from astropy import units as u

from poliastro.bodies import Earth, Mars, Sun
from poliastro.ephem import Ephem
from poliastro.frames import Planes
from poliastro.plotting import OrbitPlotter3D
from poliastro.util import time_range

EPOCH = Time("2018-02-18 12:00:00", scale="tdb")
```

```python
import plotly.io as pio
pio.renderers.default = "notebook_connected"
```

```python
roadster = Ephem.from_horizons(
    "SpaceX Roadster",
    epochs=time_range(EPOCH, end=EPOCH + 360 * u.day),
    attractor=Sun,
    plane=Planes.EARTH_ECLIPTIC,
    id_type="majorbody",
)
roadster
```

```python
from poliastro.plotting.misc import plot_solar_system
```

```python tags=["nbsphinx-thumbnail"]
frame = plot_solar_system(outer=False, epoch=EPOCH)
frame.plot_ephem(roadster, EPOCH, label="SpaceX Roadster", color="black");
```

```python
frame = OrbitPlotter3D(plane=Planes.EARTH_ECLIPTIC)

frame.plot_body_orbit(Earth, EPOCH)
frame.plot_body_orbit(Mars, EPOCH)

frame.plot_ephem(roadster, EPOCH, label="SpaceX Roadster", color="black")

frame.set_view(45 * u.deg, -120 * u.deg, 4 * u.km)
```
