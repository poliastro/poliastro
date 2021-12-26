import numpy as np
from scipy.interpolate import interp1d

__all__ = ["interp1d", "spline_interp", "sinc_interp"]


def spline_interp(y, x, u, *, kind="cubic"):
    """Interpolates y, sampled at x instants, at u instants using `scipy.interpolate.interp1d`."""
    y_u = interp1d(x, y, kind=kind)(u)
    return y_u


def sinc_interp(y, x, u):
    """Interpolates y, sampled at x instants, at u instants using sinc interpolation.

    Notes
    -----
    Taken from https://gist.github.com/endolith/1297227.
    Possibly equivalent to `scipy.signal.resample`,
    see https://mail.python.org/pipermail/scipy-user/2012-January/031255.html.
    However, quick experiments show different ringing behavior.

    """
    if len(y) != len(x):
        raise ValueError("x and s must be the same length")

    # Find the period and assume it's constant
    T = x[1] - x[0]

    sincM = np.tile(u, (len(x), 1)) - np.tile(x[:, np.newaxis], (1, len(u)))
    y_u = y @ np.sinc(sincM / T)

    return y_u
