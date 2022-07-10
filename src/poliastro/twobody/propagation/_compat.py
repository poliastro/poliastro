import sys


# https://stackoverflow.com/a/48100440
class OldPropagatorModule(sys.modules[__name__].__class__):
    def __call__(self):
        raise ImportError(
            "Propagator functions are gone, use *Propagator classes instead. "
            "This error message will be removed in a future release "
            "and you will get a 'TypeError: 'module' object is not callable'."
        )


def propagate(*args, **kwargs):
    raise ImportError(
        "The 'propagate' function is gone, use either "
        "Orbit.propagate if you want to retrieve the final state, "
        "or Orbit.to_ephem if you want to compute all the trajectory. "
        "This error message will be removed in a future release "
        "and you will get a 'ImportError: cannot import name 'propagate' "
        "from 'poliastro.twobody.propagation''."
    )
