""" Holds ground-track plotter for Earth satellites """

import plotly.graph_objects as go
from astropy import units as u
from astropy.coordinates import (
    GCRS,
    ITRS,
    CartesianDifferential,
    CartesianRepresentation,
    SphericalRepresentation,
)

from poliastro.bodies import Earth
from poliastro.earth.plotting.utils import EARTH_PALETTE
from poliastro.twobody.sampling import EpochsArray


class GroundtrackPlotter:
    """Generates two-dimensional ground-track"""

    def __init__(self, fig=None, color_palette=EARTH_PALETTE):
        """Initializes the ground-track

        Parameters
        ----------
        fig : ~plotly.graph_objects.Figure
            Figure instance for the canvas
        color_palette : dict
            A color palette for background map

        """

        # Generate custom figure if required
        if not fig:
            self.fig = go.Figure(go.Scattergeo())
        else:
            self.fig = fig

        # Default configuration is applied
        self.update_geos(
            showcoastlines=True,
            coastlinecolor="Black",
            showland=True,
            landcolor=color_palette["land_color"],
            showocean=True,
            oceancolor=color_palette["ocean_color"],
            showlakes=False,
            showrivers=False,
            lataxis={"showgrid": True, "gridcolor": "black"},
            lonaxis={"showgrid": True, "gridcolor": "black"},
        )

    def update_geos(self, **config):
        """Enables user to customize geo figure

        Parameters
        ----------
        **config : dict
            A collection of custom values for geo figure

        """
        self.fig.update_geos(config)
        return self.fig

    def update_layout(self, **config):
        """Enables user to customize figure layout

        Parameters
        ----------
        **config : dict
            A collection of custom values for figure layout

        """
        self.fig.update_layout(config)
        return self.fig

    def add_trace(self, trace):
        """Adds trace to custom figure"""

        self.fig.add_trace(trace)

    def _get_raw_coords(self, orb, t_deltas):
        """Generates raw orbit coordinates for given epochs

        Parameters
        ----------
        orb : ~poliastro.twobody.Orbit
            Orbit to be propagated
        t_deltas : ~astropy.time.DeltaTime
            Desired observation time

        Returns
        -------
        raw_xyz : numpy.ndarray
            A collection of raw cartesian position vectors
        raw_epochs : numpy.ndarray
            Associated epoch with previously raw coordinates
        """

        # Solve for raw coordinates and epochs
        ephem = orb.to_ephem(EpochsArray(orb.epoch + t_deltas))
        rr, vv = ephem.rv()
        raw_xyz = CartesianRepresentation(
            rr,
            xyz_axis=-1,
            differentials=CartesianDifferential(vv, xyz_axis=-1),
        )
        raw_epochs = ephem.epochs

        return raw_xyz, raw_epochs

    def _from_raw_to_ITRS(self, raw_xyz, raw_obstime):
        """Converts raw coordinates to ITRS ones

        Parameters
        ----------
        raw_xyz : numpy.ndarray
            A collection of rwa position coordinates
        raw_obstime : numpy.ndarray
            Associated observation time

        Returns
        -------
        itrs_xyz: ~astropy.coordinates.ITRS
            A collection of coordinates in ITRS frame

        """

        # Build GCRS and ITRS coordinates
        gcrs_xyz = GCRS(
            raw_xyz,
            obstime=raw_obstime,
            representation_type=CartesianRepresentation,
        )
        itrs_xyz = gcrs_xyz.transform_to(ITRS(obstime=raw_obstime))

        return itrs_xyz

    def _trace_groundtrack(self, orb, t_deltas, label, line_style):
        """Generates a trace for EarthSatellite's orbit grountrack

        Parameters
        ----------
        orb : ~poliastro.twobody.Orbit
            EarthSatellite's associated Orbit
        t_deltas : ~astropy.time.DeltaTime
            Collection of epochs
        label : string
            Name for the trace
        line_style : dict
            Dictionary for customizing groundtrack line trace

        Returns
        -------
        gnd_trace: ~plotly.graph_objects.Scattergeo
            Trace associated to grountrack

        """

        # Compute predicted grountrack positions
        raw_xyz, raw_obstime = self._get_raw_coords(orb, t_deltas)
        itrs_xyz = self._from_raw_to_ITRS(raw_xyz, raw_obstime)
        itrs_latlon = itrs_xyz.represent_as(SphericalRepresentation)

        # Append predicted positions to map
        gnd_trace = go.Scattergeo(
            lat=itrs_latlon.lat.to(u.deg),
            lon=itrs_latlon.lon.to(u.deg),
            mode="lines",
            name=label,
            line=line_style,
        )

        return gnd_trace

    def _trace_position(self, ss, label, marker):
        """Adds marker trace to self figure showing current position

        Parameters
        ----------
        ss : ~poliastro.twobody.Orbit
            EarthSatellite's orbit
        label : string
            Label for the orbit
        marker : dict
            Dicitonary holding plotly marker configuration

        Returns
        -------
        trace: ~plotly.graph_objects.Scattergeo
            Scattergeo trace for current position

        """

        # Check if marker available
        if not marker:
            marker = {"size": 5}

        # Solve for actual position within groundtrack
        raw_pos, raw_epoch = ss.rv()[0], ss.epoch
        itrs_pos = self._from_raw_to_ITRS(raw_pos, raw_epoch)
        itrs_latlon_pos = itrs_pos.represent_as(SphericalRepresentation)

        # Append predicted positions to map
        trace = go.Scattergeo(
            lat=itrs_latlon_pos.lat.to(u.deg),
            lon=itrs_latlon_pos.lon.to(u.deg),
            name=label,
            marker=marker,
            showlegend=False,
        )

        return trace

    def plot(self, earth_orb, t_span, label, color, line_style={}, marker={}):
        """Plots desired Earth satellite orbit for a given time span.

        Parameters
        ----------
        earth_orb : ~poliastro.earth.EarthSatellite
            Desired Earth's satellite to who's grountrack will be plotted
        t_span : ~astropy.time.TimeDelta
            A collection of epochs
        label : str
            Label for the groundtrack.
        color : string
            Desired lines and traces color
        line_style : dict
            Dictionary for customizing groundtrack line trace
        marker : dict
            Dictionary for customizing groundtrack marker trace

        Returns
        -------
        fig: ~plotly.graph_objects.Figure
            Output figure

        """

        # Retrieve basic parameters and check for proper attractor
        orb = earth_orb.orbit
        if orb.attractor != Earth:
            raise ValueError(
                f"Satellite should be orbiting Earth, not {orb.attractor}."
            )
        else:
            t_deltas = t_span - orb.epoch

        # Ensure same line and marker color unless user specifies
        for style in [line_style, marker]:
            style.setdefault("color", color)

        # Generate groundtrack trace and add it to figure
        gnd_trace = self._trace_groundtrack(orb, t_deltas, label, line_style)
        self.add_trace(gnd_trace)

        # Generate position trace and add it to figure
        pos_trace = self._trace_position(orb, label, marker)
        self.add_trace(pos_trace)

        # Return figure
        return self.fig
