import numpy as np
from numba import njit as jit


@jit
def norm(arr):
    return np.sqrt(arr @ arr)
