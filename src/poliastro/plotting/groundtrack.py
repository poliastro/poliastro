""" This script holds the function for plotting the groundtracks. """
import numpy as np
import plotly.graph_objects as go
from astropy import units as u
from astropy.coordinates import (
    GCRS,
    ITRS,
    CartesianRepresentation,
    SphericalRepresentation,
)
from astropy.time import TimeDelta

from poliastro.examples import iss
from poliastro.twobody import propagation


def groundtrack(
    orbit,
    tof=None,
    values=200,
    label="Satellite Groundtrack",
    color="blue",
    width=2,
    title="Groundtrack plot of the satellite",
    showland=True,
    showcountries=False,
    showocean=True,
    landcolor="rgb(229, 236, 246)",
    oceancolor="rgb(255, 255, 255)",
    projection="equirectangular",
    lat_grid=False,
    lat_width=0.5,
    lon_grid=False,
    lon_width=0.5,
):
    """ Plots the groundtrack of an Orbit.
        Parameters
        ----------
        orbit: poliastro.twobody.orbit.Orbit
        Orbit for making the groundtrack

        tof: float
        Time of flight in hrs

        values: int
        No. of sample points for the orbit

        label: string
        Label of the plot

        color: string of the color/rgb
        Color of the plot

        width: int
        Width of the line of the plot

        title: string
        Title of the map

        showland: True/False
        Whether to show the land or not

        showcountries: True/False
        Whether to show the countries or not

        showocean: True/False
        Whether to show the ocean or not

        landcolor: string of the color/rgb
        Color of the land

        oceancolor:  string of the color/rgb
        Color of the ocean

        projection: 'equirectangular', 'mercator', 'orthographic', 'natural earth', 'kavrayskiy7', 'miller', 'robinson', 'eckert4', 'azimuthal equal area', 'azimuthal equidistant', 'conic equal area', 'conic conformal', 'conic equidistant', 'gnomonic', 'stereographic', 'mollweide', 'hammer', 'transverse mercator', 'albers usa', 'winkel tripel', 'aitoff' and 'sinusoidal'.
        Projection type of the map

        lat_grid: True/False
        To show grid lines of the latitude or not

        lat_width: int
        Specify the width of the lat_grid

        lon_grid: True/False
        To show grid lines of the longitude or not

        lon_width: int
        Specify the width of the lon_grid
        ----------

        """

    # Calculation of time
    if tof is not None:
        time_values = TimeDelta(np.arange(1, values) * tof / values * u.h)
    else:
        time_values = TimeDelta(np.arange(1, values) * orbit.period / values)

    orbit_gcrs = propagation.propagate(orbit, time_values)

    time = {}
    time = orbit.epoch + time_values

    # Conversion of GCRS into ITRS
    orbit_itrs = GCRS(
        orbit_gcrs, obstime=time, representation_type=CartesianRepresentation
    ).transform_to(ITRS(obstime=time))

    # Converting into lat and lon
    latlon = orbit_itrs.represent_as(SphericalRepresentation)

    # Plotting the groundtrack
    fig = go.Figure()
    fig.add_trace(
        go.Scattergeo(
            lat=latlon.lat.to(u.deg),
            lon=latlon.lon.to(u.deg),
            mode="lines",
            name=label,
            showlegend=True,
            line=dict(width=width, color=color,),
        )
    )

    fig.update_layout(
        title_text=title,
        geo=dict(
            showland=showland,
            showcountries=showcountries,
            showocean=showocean,
            landcolor=landcolor,
            oceancolor=oceancolor,
            projection=dict(type=projection,),
            lonaxis=dict(showgrid=lat_grid, gridwidth=lat_width),
            lataxis=dict(showgrid=lon_grid, gridwidth=lon_width),
        ),
    )

    fig.show()


if __name__ == "__main__":
    groundtrack(iss, values=500)
