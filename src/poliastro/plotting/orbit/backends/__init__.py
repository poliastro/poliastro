"""A sub-package collecting all support orbit plotter backends."""


from poliastro.plotting.orbit.backends.matplotlib_backend import (
    OrbitPlotterBackendMatplotlib2D,
)
from poliastro.plotting.orbit.backends.plotly_backend import (
    OrbitPlotterBackendPlotly2D,
)

SUPPORTED_ORBIT_PLOTTER_BACKENDS_2D = {
    "matplotlib2D": OrbitPlotterBackendMatplotlib2D,
    "plotly2D": OrbitPlotterBackendPlotly2D,
}
"""A dictionary relating 2D orbit plotter backends and their classes."""

SUPPORTED_ORBIT_PLOTTER_BACKENDS_3D = {}
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
