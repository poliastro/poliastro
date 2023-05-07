from numba import njit as jit
import numpy as np


@jit
def norm(arr):
    return np.sqrt(arr @ arr)
