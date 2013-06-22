"""Plotting utilities.

"""

import numpy as np

from poliastro import twobody

__all__ = ['plot_orbit']


def plot_orbit(k, p, ecc, inc, omega, argp, ax, nu_vals=None, *args, **kwargs):
    """Plots orbit given orbital parameters.

    Rest of the arguments are passed to ax.plot3D.
    Based on plotorb by Richard Rieber.

    """
    # TODO: Improve for hyperbolic orbits
    # TODO: Improve number of points when speed is slow
    close = False
    if nu_vals is None:
        close = True
        nu_vals = np.linspace(0, 2 * np.pi, 200)
    N = len(nu_vals)
    r_array = np.empty((N, 3))
    for ii, nu in enumerate(nu_vals):
        r, _ = twobody.coe2rv(k, p, ecc, inc, omega, argp, nu)
        r_array[ii, :] = r
    # Close orbit
    if close:
        r_array = np.append(r_array, [r_array[0]], 0)
    return ax.plot3D(r_array[:, 0], r_array[:, 1], r_array[:, 2],
                     *args, **kwargs)