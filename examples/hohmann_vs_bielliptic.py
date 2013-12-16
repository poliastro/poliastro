import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

from poliastro import twobody
from poliastro.constants import k_Earth

ZOOM = True

v_i = 7.905  # km/s, immaterial

R = np.linspace(2, 75, num=1000)
Rstar = [15.58, 40, 60, 100, 200, np.inf]

hohmann_data = np.zeros_like(R)
bielliptic_data = np.zeros((len(R), len(Rstar)))

r_i = k_Earth / v_i ** 2
for ii, r in enumerate(R):
    r_f = r * r_i
    dva, dvb, _, _ = twobody.hohmann(k_Earth, r_i, r_f)
    hohmann_data[ii] = (abs(dva) + abs(dvb)) / v_i
    for jj, rstar in enumerate(Rstar):
        r_b = rstar * r_i
        dva, dvb, dvc, *_ = twobody.bielliptic(k_Earth, r_i, r_b, r_f)
        bielliptic_data[ii, jj] = (abs(dva) + abs(dvb) + abs(dvc)) / v_i

idx_max = np.argmax(hohmann_data)
hohmann_max = hohmann_data[idx_max]

ylims = (0.35, 0.6)

fig, ax = plt.subplots()

l, = ax.plot(R, hohmann_data, lw=2)
for jj in range(len(Rstar)):
    ax.plot(R, bielliptic_data[:, jj], color=l.get_color())
ax.vlines([11.94, R[idx_max]], *ylims, color='0.6')

if ZOOM:
    ax_zoom = zoomed_inset_axes(ax, 4, loc=4, axes_kwargs={'axisbg': '0.97'})
    ax_zoom.plot(R, hohmann_data, lw=2)
    for jj in range(len(Rstar)):
        ax_zoom.plot(R, bielliptic_data[:, jj], color=l.get_color())
    ax_zoom.vlines([11.94, R[idx_max]], *ylims, color='0.6')

    ax_zoom.set_xlim(11.0, 16.0)
    ax_zoom.set_ylim(0.52, 0.545)
    ax_zoom.set_xticklabels([], visible=False)
    ax_zoom.set_yticklabels([], visible=False)
    ax_zoom.grid(False)
    mark_inset(ax, ax_zoom, loc1=1, loc2=3, fc="none", ec='0.3')

ax.set_xlabel("R")
ax.set_ylabel("Relative change in velocity")
ax.set_ylim(*ylims)
ax.set_xlim(2, 75)
ax.set_title("Hohmann vs bielliptic transfers")

fig.show()
