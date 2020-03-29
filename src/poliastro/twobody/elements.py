import numpy as np
from astropy import units as u


def mean_motion(k, a):
    """Mean motion given body (k) and semimajor axis (a).

    """
    return np.sqrt(k / abs(a ** 3)).to(1 / u.s) * u.rad


def period(k, a):
    """Period given body (k) and semimajor axis (a).

    """
    n = mean_motion(k, a)
    return 2 * np.pi * u.rad / n
