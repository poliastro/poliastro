"""Utility hypergeometric functions.

"""
import numpy as np
from numba import njit as jit


@jit
def hyp2f1b(x):
    """Hypergeometric function 2F1(3, 1, 5/2, x), see [Battin].

    .. todo::
        Add more information about this function

    Note
    ----
    More information about hypergeometric function can be checked at
    https://en.wikipedia.org/wiki/Hypergeometric_function

    """
    if x >= 1.0:
        return np.inf
    else:
        res = 1.0
        term = 1.0
        ii = 0
        while True:
            term = term * (3 + ii) * (1 + ii) / (5 / 2 + ii) * x / (ii + 1)
            res_old = res
            res += term
            if res_old == res:
                return res
            ii += 1
