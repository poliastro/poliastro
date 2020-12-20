"""Just-in-time compiler.

Wraps numba if it is available as a module, uses an identity
decorator instead.

"""
import inspect
import warnings


def ijit(first=None, *args, **kwargs):
    """Identity JIT, returns unchanged function."""

    def _jit(f):
        return f

    if inspect.isfunction(first):
        return first
    else:
        return _jit


try:
    import numba

    jit = numba.njit
    prange = numba.prange
except ImportError:
    warnings.warn(
        "Could not import numba package. All poliastro "
        "functions will work properly but the CPU intensive "
        "algorithms will be slow. Consider installing numba to "
        "boost performance."
    )
    jit = ijit
    prange = range
