"""Just-in-time compiler.

Wraps numba if it is available as a module, uses an identity
decorator instead.

"""
import warnings
import inspect
import astropy.units as u


def ijit(first=None, *args, **kwargs):
    """Identity JIT, returns unchanged function.

    """
    def _jit(f):
        return f

    if inspect.isfunction(first):
        return first
    else:
        return _jit


try:
    import numba
    jit = numba.njit
except ImportError:
    warnings.warn("Could not import numba package. All poliastro "
                  "functions will work properly but the CPU intensive "
                  "algorithms will be slow. Consider installing numba to "
                  "boost performance.")
    jit = ijit


def accel_angles(func):
    func_fast = jit(func)

    def wrapper(angle, ecc):
        if hasattr(angle, "unit"):
            angle = angle.to(u.rad).value
        if hasattr(ecc, "unit"):
            ecc = ecc.to(u.one).value
        return func_fast(angle, ecc) * u.rad

    return wrapper
