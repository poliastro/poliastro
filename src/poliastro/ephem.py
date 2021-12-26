from enum import Enum, auto
from warnings import warn

import numpy as np
from astropy import units as u
from astropy.coordinates import (
    ICRS,
    BarycentricMeanEcliptic,
    CartesianDifferential,
    CartesianRepresentation,
    get_body_barycentric_posvel,
)
from astropy.time import Time
from astroquery.jplhorizons import Horizons
from scipy.interpolate import interp1d

from poliastro.bodies import Earth
from poliastro.frames import Planes
from poliastro.frames.util import get_frame
from poliastro.twobody.propagation import propagate
from poliastro.util import time_range
from poliastro.warnings import TimeScaleWarning

EPHEM_FORMAT = (
    "Ephemerides at {num} epochs from {start} ({start_scale}) to {end} ({end_scale})"
)


class InterpolationMethods(Enum):
    """Interpolation methods to sample the ephemerides."""

    SINC = auto()
    SPLINES = auto()


def build_ephem_interpolant(body, period, t_span, rtol=1e-5, attractor=Earth):
    """Interpolates ephemerides data

    Parameters
    ----------
    body : Body
        Source body.
    period : ~astropy.units.Quantity
        Orbital period.
    t_span : list(~astropy.units.Quantity)
        Initial and final epochs.
    rtol : float, optional
        Relative tolerance. Controls the number of sampled data points,
        defaults to 1e-5.

    Returns
    -------
    interpolant : callable
        Interpolant function that receives time increment in seconds
        since the initial epoch.

    """
    epochs = time_range(
        Time(t_span[0], format="jd", scale="tdb"),
        end=Time(t_span[1], format="jd", scale="tdb"),
        periods=int(((t_span[1] - t_span[0]) / (period * rtol)).to_value(u.one)),
    )
    ephem = Ephem.from_body(body, epochs, attractor=attractor)

    def interpolant(t):
        coords = ephem.sample(epochs[0] + (t << u.s))
        return coords.xyz[:, 0].to_value(u.km)

    return interpolant


def _interpolate_splines(epochs, reference_epochs, coordinates, *, kind="cubic"):
    xyz_unit = coordinates.xyz.unit
    d_xyz_unit = coordinates.differentials["s"].d_xyz.unit

    # TODO: Avoid building interpolant every time?
    # Does it have a performance impact?
    interpolant_xyz = interp1d(reference_epochs.jd, coordinates.xyz.value, kind=kind)
    interpolant_d_xyz = interp1d(
        reference_epochs.jd,
        coordinates.differentials["s"].d_xyz.value,
        kind=kind,
    )

    result_xyz = interpolant_xyz(epochs.jd) * xyz_unit
    result_d_xyz = interpolant_d_xyz(epochs.jd) * d_xyz_unit

    return CartesianRepresentation(
        result_xyz, differentials=CartesianDifferential(result_d_xyz)
    )


def _interpolate_sinc(epochs, reference_epochs, coordinates):
    def _interp_1d(arr):
        return _sinc_interp(arr, reference_epochs.jd, epochs.jd)

    xyz_unit = coordinates.xyz.unit
    d_xyz_unit = coordinates.differentials["s"].d_xyz.unit

    x = _interp_1d(coordinates.x.value) * xyz_unit
    y = _interp_1d(coordinates.y.value) * xyz_unit
    z = _interp_1d(coordinates.z.value) * xyz_unit

    d_x = _interp_1d(coordinates.differentials["s"].d_x.value) * d_xyz_unit
    d_y = _interp_1d(coordinates.differentials["s"].d_y.value) * d_xyz_unit
    d_z = _interp_1d(coordinates.differentials["s"].d_z.value) * d_xyz_unit

    return CartesianRepresentation(
        x, y, z, differentials=CartesianDifferential(d_x, d_y, d_z)
    )


# Taken from https://gist.github.com/endolith/1297227
def _sinc_interp(x, s, u):
    """Interpolates x, sampled at "s" instants, at "u" instants."""
    if len(x) != len(s):
        raise ValueError("x and s must be the same length")

    # Find the period and assume it's constant
    T = s[1] - s[0]

    sincM = np.tile(u, (len(s), 1)) - np.tile(s[:, np.newaxis], (1, len(u)))
    y = x @ np.sinc(sincM / T)

    return y


_INTERPOLATION_MAPPING = {
    InterpolationMethods.SINC: _interpolate_sinc,
    InterpolationMethods.SPLINES: _interpolate_splines,
}


def _get_destination_frame(attractor, plane, epochs):
    if attractor is not None:
        destination_frame = get_frame(attractor, plane, epochs)
    elif plane is Planes.EARTH_ECLIPTIC:
        destination_frame = BarycentricMeanEcliptic()
    else:
        destination_frame = None

    return destination_frame


class Ephem:
    """Time history of position and velocity of some object at particular epochs.

    Instead of creating Ephem objects directly, use the available classmethods.

    Parameters
    ----------
    coordinates : astropy.coordinates.CartesianRepresentation
        Coordinates with velocities.
    epochs : astropy.time.Time
        Epochs corresponding to the coordinates.
    plane : ~poliastro.frames.Planes
        Reference plane of the coordinates.

    """

    def __init__(self, coordinates, epochs, plane):
        if coordinates.ndim != 1 or epochs.ndim != 1:
            raise ValueError(
                f"Coordinates and epochs must have dimension 1, got {coordinates.ndim} and {epochs.ndim}"
            )

        self._epochs = epochs
        self._coordinates = coordinates
        self._plane = Planes(plane)

    def __str__(self):
        return EPHEM_FORMAT.format(
            num=len(self.epochs),
            start=self.epochs[0],
            start_scale=self.epochs[0].scale.upper(),
            end=self.epochs[-1],
            end_scale=self.epochs[-1].scale.upper(),
        )

    def __repr__(self):
        return self.__str__()

    @property
    def epochs(self):
        """Epochs at which the ephemerides was originally sampled."""
        return self._epochs

    @property
    def plane(self):
        """Reference plane of the ephemerides."""
        return self._plane

    @classmethod
    def from_body(cls, body, epochs, *, attractor=None, plane=Planes.EARTH_EQUATOR):
        """Return `Ephem` for a `SolarSystemPlanet` at certain epochs.

        Parameters
        ----------
        body: ~poliastro.bodies.SolarSystemPlanet
            Body.
        epochs: ~astropy.time.Time
            Epochs to sample the body positions.
        attractor : ~poliastro.bodies.SolarSystemPlanet, optional
            Body to use as central location,
            if not given the Solar System Barycenter will be used.
        plane : ~poliastro.frames.Planes, optional
            Fundamental plane of the frame, default to Earth Equator.

        """
        if epochs.isscalar:
            epochs = epochs.reshape(1)

        if epochs.scale != "tdb":
            epochs = epochs.tdb
            warn(
                "Input time was converted to scale='tdb' with value "
                f"{epochs.tdb.value}. Use Time(..., scale='tdb') instead.",
                TimeScaleWarning,
                stacklevel=2,
            )

        r, v = get_body_barycentric_posvel(body.name, epochs)
        coordinates = r.with_differentials(v.represent_as(CartesianDifferential))

        destination_frame = _get_destination_frame(attractor, plane, epochs)

        if destination_frame is not None:
            coordinates = (
                ICRS(coordinates)
                .transform_to(destination_frame)
                .represent_as(CartesianRepresentation, CartesianDifferential)
            )

        return cls(coordinates, epochs, plane)

    @classmethod
    def from_horizons(
        cls,
        name,
        epochs,
        *,
        attractor=None,
        plane=Planes.EARTH_EQUATOR,
        id_type=None,
    ):
        """Return `Ephem` for an object using JPLHorizons module of Astroquery.

        Parameters
        ----------
        name : str
            Name of the body to query for.
        epochs: ~astropy.time.Time
            Epochs to sample the body positions.
        attractor : ~poliastro.bodies.SolarSystemPlanet, optional
            Body to use as central location,
            if not given the Solar System Barycenter will be used.
        plane : ~poliastro.frames.Planes, optional
            Fundamental plane of the frame, default to Earth Equator.
        id_type : NoneType or str, optional
            Use "smallbody" for Asteroids and Comets and None (default) to first
            search for Planets and Satellites.

        """
        if epochs.isscalar:
            epochs = epochs.reshape(1)

        refplanes_dict = {
            Planes.EARTH_EQUATOR: "earth",
            Planes.EARTH_ECLIPTIC: "ecliptic",
        }
        refplane = refplanes_dict[plane]

        if attractor is not None:
            bodies_dict = {
                "sun": 10,
                "mercury": 199,
                "venus": 299,
                "earth": 399,
                "mars": 499,
                "jupiter": 599,
                "saturn": 699,
                "uranus": 799,
                "neptune": 899,
            }
            location = f"500@{bodies_dict[attractor.name.lower()]}"
        else:
            location = "@ssb"

        obj = Horizons(
            id=name, location=location, epochs=epochs.jd, id_type=id_type
        ).vectors(refplane=refplane)

        x = obj["x"]
        y = obj["y"]
        z = obj["z"]
        d_x = obj["vx"]
        d_y = obj["vy"]
        d_z = obj["vz"]

        coordinates = CartesianRepresentation(
            x, y, z, differentials=CartesianDifferential(d_x, d_y, d_z)
        )
        return cls(coordinates, epochs, plane)

    @classmethod
    def from_orbit(
        cls,
        orbit,
        epochs,
        plane=Planes.EARTH_EQUATOR,
    ):
        """Return `Ephem` from an `Orbit` at certain epochs.

        Parameters
        ----------
        orbit: ~poliastro.twobody.orbit.Orbit
            Orbit.
        epochs: ~astropy.time.Time
            Epochs to sample the orbit positions.
        plane: ~poliastro.frames.Planes, optional
            Fundamental plane of the frame, default to Earth Equator.

        """
        if epochs.isscalar:
            epochs = epochs.reshape(1)

        time_of_flights = epochs - orbit.epoch
        coordinates = propagate(orbit, time_of_flights)

        return cls(coordinates, epochs, plane)

    def sample(self, epochs=None, *, method=InterpolationMethods.SPLINES, **kwargs):
        """Returns coordinates at specified epochs.

        Parameters
        ----------
        epochs: ~astropy.time.Time, optional
            Epochs to sample the ephemerides,
            if not given the original one from the object will be used.
        method : ~poliastro.ephem.InterpolationMethods, optional
            Interpolation method to use for epochs outside of the original ones,
            default to splines.

        Returns
        -------
        CartesianRepresentation
            Sampled coordinates with velocities.

        """
        if epochs is None or epochs.isscalar and (epochs == self.epochs).all():
            return self._coordinates

        # TODO: Proper type annotation
        coordinates = _INTERPOLATION_MAPPING[method](  # type: ignore
            epochs.reshape(-1), self.epochs, self._coordinates, **kwargs
        )

        return coordinates

    def rv(self, epochs=None, **kwargs):
        """Position and velocity vectors at given epochs.

        Parameters
        ----------
        epochs : ~astropy.time.Time, optional
            Epochs to sample the ephemerides, default to now.

        """
        coordinates = self.sample(epochs, **kwargs)

        r = coordinates.get_xyz(xyz_axis=1)
        v = coordinates.differentials["s"].get_d_xyz(xyz_axis=1)

        if epochs is not None and epochs.isscalar:
            return r[0], v[0]
        else:
            return r, v
