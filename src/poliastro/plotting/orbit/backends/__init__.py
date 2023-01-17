"""A sub-package collecting all support orbit plotter backends."""


from poliastro.plotting.orbit.backends.matplotlib import Matplotlib2D
from poliastro.plotting.orbit.backends.plotly import Plotly2D, Plotly3D

DEFAULT_ORBIT_PLOTTER_BACKENDS_2D = {
    "matplotlib2D": Matplotlib2D,
    "plotly2D": Plotly2D,
}
"""A dictionary relating 2D orbit plotter backends and their classes."""

DEFAULT_ORBIT_PLOTTER_BACKENDS_3D = {
    "plotly3D": Plotly3D,
}
"""A dictionary relating 3D orbit plotter backends and their classes."""


DEFAULT_ORBIT_PLOTTER_BACKENDS = {
    **DEFAULT_ORBIT_PLOTTER_BACKENDS_2D,
    **DEFAULT_ORBIT_PLOTTER_BACKENDS_3D,
}
"""A dictionary relating orbit plotter backends and their classes."""

DEFAULT_ORBIT_PLOTTER_MATPLOTLIB_BACKENDS = {
    "matplotlib2D": Matplotlib2D,
}
"""A dictionary relating orbit plotter backends for matplotlib and their classes."""

DEFAULT_ORBIT_PLOTTER_PLOTLY_BACKENDS = {
    "plotly2D": Plotly2D,
    "plotly3D": Plotly3D,
}
"""A dictionary relating orbit plotter backends for plotly and their classes."""

__all__ = [
    "DEFAULT_ORBIT_PLOTTER_BACKENDS_2D",
    "DEFAULT_ORBIT_PLOTTER_BACKENDS_3D",
    "DEFAULT_ORBIT_PLOTTER_BACKENDS",
    "DEFAULT_ORBIT_PLOTTER_MATPLOTLIB_BACKENDS",
    "DEFAULT_ORBIT_PLOTTER_PLOTLY_BACKENDS",
]
