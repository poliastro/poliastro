import numpy as np
from astropy import units as u
from matplotlib import animation, patches as mpl_patches, pyplot as plt
from matplotlib.colors import to_rgba

from ._base import BaseOrbitPlotter, Mixin2D


class AnimatedOrbitPlotter(BaseOrbitPlotter, Mixin2D):
    """AnimatedOrbitPlotter class.

    This class holds the perifocal plane of the first
    :py:class:`~poliastro.twobody.orbit.Orbit` plotted in it using
    :py:meth:`plot`, so all following
    plots will be projected on that plane. Alternatively, you can call
    :py:meth:`set_frame` to set the frame before plotting.

    """

    def __init__(self, ax=None, fig=None, num_points=150, dark=False, *, plane=None):
        """Constructor.

        Parameters
        ----------
        ax : ~matplotlib.axes.Axes
            Axes in which to plot. If not given, new ones will be created.
        num_points : int, optional
            Number of points to use in plots, default to 150.
        dark : bool, optional
            If set as True, plots the orbit in Dark mode.

        """
        super().__init__(num_points=num_points, plane=plane)

        self._ax = ax
        self.fig = fig
        if not self._ax:
            if dark:
                with plt.style.context("dark_background"):
                    self.fig, self._ax = plt.subplots(figsize=(6, 6))
            else:
                self.fig, self._ax = plt.subplots(figsize=(6, 6))

        self._frame = None

    def _clear_attractor(self):
        for attractor in self._ax.findobj(match=mpl_patches.Circle):
            attractor.remove()

    def _get_colors(self, color, trail):
        if color is None:
            # HACK: https://stackoverflow.com/a/13831816/554319
            color = next(self._ax._get_lines.prop_cycler)["color"]

        if trail:
            colors = [color, to_rgba(color, 0)]
        else:
            colors = [color]
        return colors

    def _draw_point(self, radius, color, name, center=None):
        x_center, y_center = self._project(
            center[None]
        )  # Indexing trick to add one extra dimension

        (l,) = self._ax.plot(
            x_center.to(u.km).value, y_center.to(u.km).value, "o", mew=0, color=color
        )

        return l

    def _draw_sphere(self, radius, color, name, center=[0, 0, 0] * u.km):
        x_center, y_center = self._project(
            center[None]
        )  # Indexing trick to add one extra dimension

        self._ax.add_patch(
            mpl_patches.Circle(
                (x_center.to(u.km).value, y_center.to(u.km).value),
                radius.to(u.km).value,
                lw=0,
                color=color,
            )
        )

    def _plot_coordinates_anim(self, coordinates, label, colors, dashed):
        if dashed:
            linestyle = "dashed"

        rr = coordinates.xyz.transpose()
        x, y = self._project(rr)

        def animate(i):

            Element.set_data(x.to(u.km).value[i], y.to(u.km).value[i])
            return (Element,)

        self._ax.plot(
            x.to(u.km).value, y.to(u.km).value, linestyle=linestyle, color=colors[0]
        )
        (Element,) = self._ax.plot(
            x.to(u.km).value[0], y.to(u.km).value[0], "o", mew=0, color=colors[0]
        )
        Animation = animation.FuncAnimation(
            self.fig,
            animate,
            frames=np.arange(0, len(x), 1),
            interval=30,
            blit=True,
            repeat=True,
        )
        self._ax.legend(
            [label],
            loc="upper left",
            bbox_to_anchor=(1.05, 1.015),
            title="Names and initial epoch",
            numpoints=1,
        )

        self._ax.set_xlabel("$x$ (km)")
        self._ax.set_ylabel("$y$ (km)")
        self._ax.set_aspect(1)
        return Animation

    def anim(self, orbit, *, label=None, color=None, trail=False):
        """animates state and osculating orbit in their plane.

        Parameters
        ----------
        orbit : ~poliastro.twobody.orbit.Orbit
            Orbit to plot.
        label : string, optional
            Label of the orbit.
        color : string, optional
            Color of the line and the position.
        trail : bool, optional
            Fade the orbit trail, default to False.

        """
        if not self._frame:
            self.set_orbit_frame(orbit)

        _animation = self._anim(orbit, label=label, color=color, trail=trail)
        return _animation
