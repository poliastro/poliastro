"""A sub-package collecting all support orbit plotter backends."""


from poliastro.plotting.orbit.backends.matplotlib import (
    OrbitPlotterBackendMatplotlib2D,
)
from poliastro.plotting.orbit.backends.plotly import (
    OrbitPlotterBackendPlotly2D,
    OrbitPlotterBackendPlotly3D,
)

SUPPORTED_ORBIT_PLOTTER_BACKENDS_2D = {
    "matplotlib2D": OrbitPlotterBackendMatplotlib2D,
    "plotly2D": OrbitPlotterBackendPlotly2D,
}
"""A dictionary relating 2D orbit plotter backends and their classes."""

SUPPORTED_ORBIT_PLOTTER_BACKENDS_3D = {
    "plotly3D": OrbitPlotterBackendPlotly3D,
}
"""A dictionary relating 3D orbit plotter backends and their classes."""


SUPPORTED_ORBIT_PLOTTER_BACKENDS = {
    **SUPPORTED_ORBIT_PLOTTER_BACKENDS_2D,
    **SUPPORTED_ORBIT_PLOTTER_BACKENDS_3D,
}
"""A dictionary relating orbit plotter backends and their classes."""


__all__ = [
    "SUPPORTED_ORBIT_PLOTTER_BACKENDS_2D",
    "SUPPORTED_ORBIT_PLOTTER_BACKENDS_3D",
    "SUPPORTED_ORBIT_PLOTTER_BACKENDS",
]
