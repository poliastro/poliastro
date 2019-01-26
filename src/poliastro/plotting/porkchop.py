"""
This is the script for porkchop plotting
"""

from poliastro.bodies import Sun
from poliastro.util import norm
from poliastro.iod import lambert
from poliastro.util import time_range


from astropy import units as u
from astropy.time import Time
from astropy import coordinates as coord

import matplotlib.pyplot as plt

import numpy as np


def _targetting(departure_body, target_body, t_launch, t_arrival):
    """
    This function returns the increment in departure and arrival velocities.
    """

    # Compute departure and arrival positions
    rr_dpt_body, vv_dpt_body = coord.get_body_barycentric_posvel(
        departure_body.name, t_launch
    )
    rr_arr_body, vv_arr_body = coord.get_body_barycentric_posvel(
        target_body.name, t_arrival
    )

    # Compute time of flight
    tof = t_arrival - t_launch

    if tof <= 0:
        return None, None, None

    try:
        (v_dpt, v_arr), = lambert(Sun.k, rr_dpt_body.xyz, rr_arr_body.xyz, tof)

        # Compute all the output variables
        dv_dpt = norm(v_dpt - vv_dpt_body.xyz)
        dv_arr = norm(v_arr - vv_arr_body.xyz)
        c3_launch = dv_dpt.value ** 2
        c3_arrival = dv_arr.value ** 2

        return c3_launch, c3_arrival, tof.jd

    except AssertionError:
        return None, None, None


# numpy.vectorize is amazing
targetting_vec = np.vectorize(
    _targetting,
    otypes=[np.ndarray, np.ndarray, np.ndarray],
    excluded=["departure_body", "target_body"],
)


def porkchop(departure_body, target_body, launch_span, arrival_span):
    """Plots porkchop between two bodies.

    Parameters
    ----------
    departure_body: poliastro.bodies.Body
        Body from which departure is done
    target_body: poliastro.bodies.Body
        Body for targetting
    launch_span: astropy.time.Time
        Time span for launch
    arrival_span: astropy.time.Time
        Time span for arrival

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
    """

    c3_launch, c3_arrival, tof = targetting_vec(
        departure_body,
        target_body,
        launch_span[np.newaxis, :],
        arrival_span[:, np.newaxis],
    )

    """
    Algorithm works: 'for each launch get all arrivals'.
    Contourf works: 'for each Y -> all X'.
    We need to transpose the arrays.
    """

    fig, ax = plt.subplots(figsize=(15, 15))
    c3_levels = np.linspace(0, 45, 30)
    t_levels = np.linspace(100, 500, 5)

    c = plt.contourf(
        [D.to_datetime() for D in launch_span],
        [A.to_datetime() for A in arrival_span],
        c3_launch,
        c3_levels,
    )

    l = plt.contour(
        [D.to_datetime() for D in launch_span],
        [A.to_datetime() for A in arrival_span],
        c3_launch,
        c3_levels,
        colors="black",
        linestyles="solid",
    )

    t = plt.contour(
        [D.to_datetime() for D in launch_span],
        [A.to_datetime() for A in arrival_span],
        tof,
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
            departure_body.name, target_body.name, launch_span[0].datetime.year
        ),
        fontsize=14,
        fontweight="bold",
    )

    plt.xlabel("Launch date", fontsize=10, fontweight="bold")
    plt.ylabel("Arrival date", fontsize=10, fontweight="bold")

    return c3_launch, c3_arrival, tof
