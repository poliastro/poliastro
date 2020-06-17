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

# Plotting in 3D

```python
import numpy as np

from poliastro.constants import J2000
from poliastro.examples import *
from poliastro.plotting import *
```

```python
import plotly.io as pio
pio.renderers.default = "notebook_connected"
```

```python
churi.plot(interactive=True, use_3d=True)
```

```python
frame = OrbitPlotter3D()

frame.plot(churi)
frame.plot_body_orbit(Earth, J2000)
```

```python
frame = OrbitPlotter3D()

frame.plot(molniya)
frame.plot(iss)
```

```python
eros = Orbit.from_sbdb("eros")

frame = OrbitPlotter3D()

frame.plot_body_orbit(Earth, J2000)
frame.plot(eros, label="eros")
```

```python
from poliastro.ephem import Ephem
from poliastro.util import time_range
```

```python
date_launch = time.Time("2011-11-26 15:02", scale="utc").tdb
date_arrival = time.Time("2012-08-06 05:17", scale="utc").tdb

earth = Ephem.from_body(
    Earth, time_range(date_launch, end=date_arrival, periods=50)
)
```

```python
frame = OrbitPlotter3D()
frame.set_attractor(Sun)

frame.plot_body_orbit(Earth, J2000, label=Earth)
frame.plot_ephem(earth, label=Earth)
```

```python
frame = OrbitPlotter3D()

frame.plot(eros, label="eros")
frame.plot_trajectory(earth.sample(), label=Earth)
```
