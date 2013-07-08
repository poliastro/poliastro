"""Plotting utilities.

"""

import numpy as np

from . import twobody

__all__ = ['plot_orbit']


def plot_orbit(ax, k, a, ecc, inc, omega, argp, nu_lims=(0, 2 * np.pi), N=200,
               *args, **kwargs):
    """Plots orbit given orbital parameters.

    The function always assumes direct motion: if nu_lims[1] < nu_lims[0],
    sums one revolution to the final position.

    The rest of the arguments are passed to ax.plot3D.

    Parameters
    ----------
    ax : Axes3D
        3D axes to plot the orbit.
    k : float
        Standard gravitational parameter (km^3 / s^2).
    a : float
        Semi-major axis (km).
    ecc : float
        Eccentricity.
    inc : float
        Inclination (rad).
    omega : float
        Longitude of ascending node (rad).
    argp : float
        Argument of perigee (rad).
    nu_lims : tuple (optional)
        Minimum and maximum values of the true anomaly. Default (0, 2 * pi)
        (plot complete orbit).
    N : int (optional)
        Number of points to plot, default to 200.

    Notes
    -----
    Based on plotorb by Richard Rieber.

    """
    # TODO: Adaptative number of points?
    # TODO: Fix for hyperbolic orbits
    N = 200
    close = False
    if nu_lims is (0, 2 * np.pi):
        close = True
    else:
        while nu_lims[1] < nu_lims[0]:
            nu_lims = (nu_lims[0], nu_lims[1] + 2 * np.pi)
    nu_vals = np.linspace(nu_lims[0], nu_lims[1], N)
    r_array = np.empty((N, 3))
    for ii, nu in enumerate(nu_vals):
        r, _ = twobody.coe2rv(k, a, ecc, inc, omega, argp, nu)
        r_array[ii, :] = r
    if close:
        r_array = np.append(r_array, [r_array[0]], 0)
    return ax.plot3D(r_array[:, 0], r_array[:, 1], r_array[:, 2],
                     *args, **kwargs)
