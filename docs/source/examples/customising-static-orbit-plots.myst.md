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

# Customising static orbit plots

+++

The default styling for plots works pretty well however sometimes you may need to change things. The following will show you how to change the style of your plots and have different types of lines and dots.

This is the default plot we will start with:

```{code-cell} ipython3
from astropy.time import Time

from matplotlib import pyplot as plt

from poliastro.bodies import Earth, Mars, Jupiter, Sun
from poliastro.frames import Planes
from poliastro.plotting import OrbitPlotter
from poliastro.twobody import Orbit
```

```{code-cell} ipython3
epoch = Time("2018-08-17 12:05:50", scale="tdb")

plotter = OrbitPlotter(plane=Planes.EARTH_ECLIPTIC)
plotter.plot_body_orbit(Earth, epoch, label="Earth")
plotter.plot_body_orbit(Mars, epoch, label="Mars")
plotter.plot_body_orbit(Jupiter, epoch, label="Jupiter")
```

```{code-cell} ipython3
:tags: [nbsphinx-thumbnail]

epoch = Time("2018-08-17 12:05:50", scale="tdb")

plotter = OrbitPlotter(plane=Planes.EARTH_ECLIPTIC)
earth_plots_traj, earth_plots_pos = plotter.plot_body_orbit(
    Earth, epoch, label=Earth
)

earth_plots_traj.set_linestyle("-")  # solid line
earth_plots_traj.set_linewidth(0.5)
earth_plots_pos.set_marker("H")  # Hexagon
earth_plots_pos.set_markersize(15)

mars_plots = plotter.plot_body_orbit(Mars, epoch, label=Mars)
jupiter_plots = plotter.plot_body_orbit(Jupiter, epoch, label=Jupiter)
```

Here we get hold of the lines list from the `OrbitPlotter.plot` method this is a list of lines. The first is the orbit line. The second is the current position marker. With the matplotlib lines objects we can start changing the style. First we make the line solid but thin. Then we change the current position marker to a large hexagon.

More details of the style options for the markers can be found here: https://matplotlib.org/2.0.2/api/markers_api.html#module-matplotlib.markers.
More details of the style options on lines can be found here: https://matplotlib.org/2.0.2/api/lines_api.html. However make sure that you use the set methods rather than just changing the attributes as the methods will force a re-draw of the plot.

Next we will make some changes to the other two orbits:

```{code-cell} ipython3
epoch = Time("2018-08-17 12:05:50", scale="tdb")

plotter = OrbitPlotter()

earth_plots_t, earth_plots_p = plotter.plot_body_orbit(
    Earth, epoch, label=Earth
)
earth_plots_t.set_linestyle("-")  # solid line
earth_plots_t.set_linewidth(0.5)
earth_plots_p.set_marker("H")  # Hexagon
earth_plots_p.set_markersize(15)

mars_plots_t, mars_plots_p = plotter.plot_body_orbit(Mars, epoch, label=Mars)
mars_plots_t.set_dashes([0, 1, 0, 1, 1, 0])
mars_plots_t.set_linewidth(2)
mars_plots_p.set_marker("D")  # Diamond
mars_plots_p.set_markersize(15)
mars_plots_p.set_fillstyle("none")
# make sure this is set if you use fillstyle 'none'
mars_plots_p.set_markeredgewidth(1)

jupiter_plots_t, jupiter_plots_p = plotter.plot_body_orbit(
    Jupiter, epoch, label=Jupiter
)
jupiter_plots_t.set_linestyle("")  # No line
jupiter_plots_p.set_marker("*")  # star
jupiter_plots_p.set_markersize(15)
```

You can also change the style of the plot using the matplotlib axis which can be aquired from the `OrbitPlotter()`.

See the following example that creates a grid, adds a title, and makes the background transparent. To make the changes clearer it goes back to the inital example:

```{code-cell} ipython3
epoch = Time("2018-08-17 12:05:50", scale="tdb")

fig, ax = plt.subplots()

ax.grid(True)
ax.set_title("Earth, Mars, and Jupiter")
ax.set_facecolor("None")

plotter = OrbitPlotter(scene=ax)

plotter.plot_body_orbit(Earth, epoch, label=Earth)
plotter.plot_body_orbit(Mars, epoch, label=Mars)
plotter.plot_body_orbit(Jupiter, epoch, label=Jupiter)
```
