"""
This is the script for porkchop plotting
"""

from poliastro.bodies import Sun
from poliastro.util import norm
from poliastro.iod import lambert

from astropy import units as u
from astropy.time import Time
from astropy import coordinates as coord

import matplotlib.pyplot as plt

import numpy as np


def lambert_porkchop(dpt_body, arr_body, dpt_t, arr_t):
    """
    This function returns the increment in departure and arrival velocities.
    """

    # Compute departure and arrival positions
    rr_dpt_body, vv_dpt_body = coord.get_body_barycentric_posvel(dpt_body.name, dpt_t)
    rr_arr_body, vv_arr_body = coord.get_body_barycentric_posvel(arr_body.name, arr_t)

    # Compute time of flight
    tof = arr_t - dpt_t

    if tof <= 0:
        return None, None, None, None, None

    try:
        (v_dpt, v_arr), = lambert(Sun.k, rr_dpt_body.xyz, rr_arr_body.xyz, tof)

        # Compute all the output variables
        dv_dpt = norm(v_dpt - vv_dpt_body.xyz)
        dv_arr = norm(v_arr - vv_arr_body.xyz)
        c3_launch = dv_dpt.value ** 2
        c3_arrival = dv_arr.value ** 2

        return dv_dpt.value, dv_arr.value, c3_launch, c3_arrival, tof.jd

    except AssertionError:
        return None, None, None, None, None


# numpy.vectorize is amazing
lambert_porkchop_vectorized = np.vectorize(
    lambert_porkchop,
    otypes=[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    excluded=["dpt_body", "arr_body"],
)


def porkchop(body_dpt, body_arr, dpt_start, dpt_end, arr_start, arr_end, N=50):
    """Plots porkchop between two bodies.

    Parameters
    ----------
    body_dpt: poliastro.bodies.Body
        Body for launch
    body_arr: poliastro.bodies.Body
        Body for arrival
    dpt_start: str
        Porkchop launch date starts in this value
    dpt_end: str
        Porkchop launch date ends in this value
    arr_start: str
        Porkchop arrival date starts in this value
    arr_end: str
        Porkchop arrival date ends in this value
    
    Returns
    -------
    dpt: np.array
        Departure time span
    arr: np.array
        Arrival time span
    deltav_dpt: np.ndarray
        Departure velocity needed for each time of flight
    deltav_arr: np.ndarray
        Arrival velocity needed for each time of flight
    c3_dpt: np.ndarray
        Characteristic launch energy
    c3_arr: np.ndarray
        Characteristic arrival energy

    Example
    -------
    # Time requirements YYYY-MM-DD
    # Data is from porkchop pag. 180
    >>> from poliastro.plotting.porkchop import porkchop
    >>> from poliastro.bodies import Earth, Mars
    >>> import matplotlib.pyplot as plt
    >>> departure_start = "2005-04-30"
    >>> departure_end   = "2005-10-07"
    >>> arrival_start = "2005-11-16"
    >>> arrival_end   = "2006-12-21"
    >>> dpt, arr, dv_dpt, dv_arr, c3dpt, c3arr = porkchop(Earth, Mars, departure_start, departure_end, arrival_start, arrival_end)
    >>> plt.show()
    """

    # Computing time spans fot departure and arrival
    dpt = [
        Time(d, format="jd")
        for d in np.linspace(Time(dpt_start).jd, Time(dpt_end).jd, N + 1)
    ]
    arr = [
        Time(d, format="jd")
        for d in np.linspace(Time(arr_start).jd, Time(arr_end).jd, N + 1)
    ]

    # Prellocate in memory the arrays
    deltav_dpt = np.zeros((len(dpt), len(arr)))
    deltav_arr = np.zeros((len(dpt), len(arr)))
    c3_dpt = np.zeros((len(dpt), len(arr)))
    c3_arr = np.zeros((len(dpt), len(arr)))
    iso_tof = np.zeros((len(dpt), len(arr)))
    idx = 0

    for d in dpt:

        dv_dpt, dv_arr, c3_d, c3_a, t_flight = lambert_porkchop_vectorized(
            body_dpt, body_arr, d, arr
        )

        deltav_dpt[idx] = dv_dpt
        deltav_arr[idx] = dv_arr
        c3_dpt[idx] = c3_d
        c3_arr[idx] = c3_a
        iso_tof[idx] = t_flight
        idx += 1

    """
    Algorithm works: 'for each launch get all arrivals'.
    Contourf works: 'for each Y -> all X'.
    We need to transpose the arrays.
    """

    fig, ax = plt.subplots(figsize=(15, 15))
    c3_levels = np.linspace(0, 45, 30)
    t_levels = np.linspace(100, 500, 5)

    c = plt.contourf(
        [D.to_datetime() for D in dpt],
        [A.to_datetime() for A in arr],
        np.transpose(c3_dpt),
        c3_levels,
    )

    l = plt.contour(
        [D.to_datetime() for D in dpt],
        [A.to_datetime() for A in arr],
        np.transpose(c3_dpt),
        c3_levels,
        colors="black",
        linestyles="solid",
    )

    t = plt.contour(
        [D.to_datetime() for D in dpt],
        [A.to_datetime() for A in arr],
        np.transpose(iso_tof),
        t_levels,
        colors="red",
        linestyles="dashed",
        linewidths=3.5,
    )

    cbar = plt.colorbar(c)
    cbar.set_label("$km^2/s^2$")
    plt.clabel(l, inline=1, fmt="%1.1f", colors="k", fontsize=10)
    plt.clabel(t, inline=1, fmt="%1.1f", colors="r", fontsize=14)

    plt.grid()
    fig.autofmt_xdate()

    plt.title(
        "{} - {} for year {}, C3 Launch, TFL".format(
            body_dpt.name, body_arr.name, dpt[0].datetime.year
        ),
        fontsize=14,
        fontweight="bold",
    )

    plt.xlabel("Launch date", fontsize=10, fontweight="bold")
    plt.ylabel("Arrival date", fontsize=10, fontweight="bold")
    plt.show()

    return dpt, arr, deltav_dpt, deltav_arr, c3_dpt, c3_arr
