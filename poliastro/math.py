# coding: utf-8
import numba


@numba.njit('f8(f8[:], f8[:])')
def dot(u, v):
    """Returns the dot product of two vectors.

    """
    dp = 0.0
    for ii in range(u.shape[0]):
        dp += u[ii] * v[ii]
    return dp
