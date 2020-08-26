"""Robust smoothing functions.

Translated from Garcia, Damien. 2010. "Robust Smoothing of Gridded Data in One and Higher Dimensions with Missing Values."
Computational Statistics & Data Analysis 54 (4): 1167–78. https://doi.org/10.1016/j.csda.2009.09.020.

"""
import numpy as np
from numpy.linalg import norm
from numpy.matlib import repmat
from scipy.fftpack import dctn, idctn
from scipy.optimize import fminbound


def dct2(x):
    """2D discrete cosine transform."""
    return dctn(x, type=2, norm="ortho")


def idct2(x):
    """2D inverse discrete cosine transform."""
    return idctn(x, type=2, norm="ortho")


def rsmooth(y):
    if y.ndim < 2:
        y = np.atleast_2d(y)
        one_dim = True
    else:
        one_dim = False

    n1, n2 = y.shape
    n = n1 * n2  # noqa: F841
    N = (np.array([n1, n2]) != 1).sum()
    Lambda = (
        repmat(-2 + 2 * np.cos(np.arange(0, n2) * np.pi / n2), n1, 1)
        + (-2 + 2 * np.cos(np.arange(0, n1) * np.pi / n1))[:, None]
    )
    W = np.ones((n1, n2))
    z = zz = y

    def GCVscore(p):
        """Generalized cross-validation score."""
        # This makes the code more similar to the original
        # and avoids recomputing z after optimizing the penalty
        nonlocal z
        n = y.size
        s = 10 ** p  # Penalty term
        Gamma = 1 / (1 + s * Lambda ** 2)  # See equation 6
        z = idct2(Gamma * DCTy)
        RSS = norm(np.sqrt(W) * (y - z)) ** 2  # Residual sum-of-squares
        TrH = np.sum(Gamma)  # Trace of "hat matrix"
        GCVs = RSS / n / (1 - TrH / n) ** 2
        return GCVs

    for k in range(1, 7):
        tol = np.inf
        while tol > 1e-5:
            DCTy = dct2(W * (y - zz) + zz)
            p = fminbound(GCVscore, -15, 38)
            tol = norm(zz - z) / norm(z)
            zz = z
        s = 10 ** p
        tmp = np.sqrt(1 + 16 * s)
        h = (np.sqrt(1 + tmp) / np.sqrt(2) / tmp) ** N
        W = bisquare(y - z, h)

    if one_dim:
        z = z[0]

    return z


def bisquare(r, h):
    """Bisquare weight function.

    Notes
    -----
    Heiberger, Richard M., and Richard A. Becker. 1992. "Design of an S Function for Robust Regression
    Using Iteratively Reweighted Least Squares." Journal of Computational and Graphical Statistics 1 (3): 181–96.
    https://doi.org/10.1080/10618600.1992.10474580.

    """
    c = 4.685  # Tuning constant for a given distribution
    MAD = np.median(np.abs(r - np.median(r)))  # Median absolute deviation
    u = np.abs(r / (1.4826 * MAD) / np.sqrt(1 - h))  # Studentized residual
    W = (1 - (u / c) ** 2) ** 2 * (
        (u / c) < 1
    )  # Trick for stepwise (1 - ...) ** 2 or 0
    return W
