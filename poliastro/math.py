# coding: utf-8
import numba


def factorial(nn):
    """Return nn factorial.

    Raises ValueError if nn is not integer, negative or results in
    64-bit integer overflow.

    """
    if not isinstance(nn, int) and not nn.is_integer():
        raise ValueError("factorial() only accepts integer values")
    elif nn < 0:
        raise ValueError("factorial() not defined for negative values")
    elif nn > 20:
        raise ValueError("factorial() too big to be stored in int-64")

    return _factorial(nn)

@numba.jit('i8(i8)', nopython=True)
def _factorial(nn):
    res = 1
    for ii in range(2, nn + 1):
        res *= ii
    return res


@numba.njit('f8(f8[:], f8[:])')
def dot(u, v):
    """Returns the dot product of two vectors.

    """
    dp = 0.0
    for ii in range(u.shape[0]):
        dp += u[ii] * v[ii]
    return dp
