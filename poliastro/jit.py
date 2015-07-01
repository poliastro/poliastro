# coding: utf-8
"""Just-in-time compiler.

Wraps numba if it is available as a module, uses an identity
decorator instead.

"""
import warnings
import inspect


def ijit(first=None, *args, **kwargs):
    """Identity JIT, returns unchanged function.

    """
    def _jit(f):
        return f

    if inspect.isfunction(first):
        return first
    else:
        return _jit


def select_jit():
    """Selects appropriate JIT function.

    Returns numba.njit (nopython JIT) if available, else returns an identity
    decorator.

    """
    try:
        import numba
        jit = numba.njit
    except ImportError:
        warnings.warn("Could not import numba package. All poliastro "
                      "functions will work properly but the CPU intensive "
                      "algorithms will be slow. Consider installing numba to "
                      "boost performance")
        jit = ijit

    return jit

jit = select_jit()
